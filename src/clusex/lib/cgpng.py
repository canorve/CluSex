#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import warnings
from pathlib import Path
import tkinter as tk
import tkinter.font as tkFont
from tkinter import messagebox as tkMessageBox
from PIL import Image, ImageTk
import pandas as pd

warnings.filterwarnings("ignore")


def natural_key(s: str):
    """Human sorting: img2, img10 -> 2, 10."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]


def cgPNG(img_folder, tsv_out, include_header=False):
    img_folder = Path(img_folder)
    # Accept common formats
    paths = [p for p in img_folder.iterdir() if p.suffix.lower() in {".png", ".jpg", ".jpeg"}]
    paths.sort(key=lambda p: natural_key(p.name))

    if not paths:
        print("No images found in path")
        return

    root = tk.Tk()
    GUI(root, paths, Path(tsv_out), include_header=include_header)
    root.mainloop()


class GUI:
    def __init__(self, master, png_paths, pngfile_out, include_header=False, extensiones=(".png", ".jpg", ".jpeg")):
        self.master = master
        self.include_header = include_header
        self.extensiones = tuple(e.lower() for e in extensiones)

        # Store as strings to avoid Path/string method mismatches
        self.all_paths = [str(p) for p in png_paths if p.suffix.lower() in self.extensiones]
        self.n = len(self.all_paths)
        self.pngfile_out = Path(pngfile_out)

        # State
        self.index = 0
        self.map_cls = {}  # {path_str: class}
        self._current_img_pil = None
        self._photo = None

        # Resume if file exists
        if self.pngfile_out.exists():
            if tkMessageBox.askyesno("Continue", f"{self.pngfile_out} file already exists. Do you want to continue where you left off?"):
                df_prev = pd.read_csv(self.pngfile_out, sep=',', header=None, names=['path', 'class'], dtype=str)
                df_prev = df_prev.dropna(subset=['path', 'class'])
                self.map_cls = dict(zip(df_prev['path'], df_prev['class']))
                # Set index to first unclassified
                for i, p in enumerate(self.all_paths):
                    if p not in self.map_cls:
                        self.index = i
                        break
                else:
                    self.index = self.n

        # --- Window and grid skeleton ---
        master.title("Galaxy classifier")
        width, height = 640, 640
        sw, sh = master.winfo_screenwidth(), master.winfo_screenheight()
        x, y = int((sw - width) / 2), int((sh - height) / 2)
        master.geometry(f"{width}x{height}+{x}+{y}")
        master.resizable(True, True)
        master.minsize(560, 560)
        master.configure(bg="#1e1e1e")

        # Grid weights: canvas expands, others natural
        master.grid_rowconfigure(0, weight=0)   # header
        master.grid_rowconfigure(1, weight=1)   # canvas (grows)
        master.grid_rowconfigure(2, weight=0)   # class buttons
        master.grid_rowconfigure(3, weight=0)   # control buttons
        master.grid_columnconfigure(0, weight=1)

        self.font = tkFont.Font(family='Helvetica', size=11, weight="bold")

        # --- Header ---
        header = tk.Frame(master, bg="#1e1e1e")
        header.grid(row=0, column=0, sticky="ew", padx=12, pady=(8, 4))
        header.grid_columnconfigure(0, weight=1)
        header.grid_columnconfigure(1, weight=1)

        self.ImgTitle = tk.Label(header, font=self.font, fg="#cfcfcf", bg="#1e1e1e", anchor="w")
        self.ImgTitle.grid(row=0, column=0, sticky="w")
        self.lbl_progress = tk.Label(header, font=self.font, fg="#cfcfcf", bg="#1e1e1e", anchor="e")
        self.lbl_progress.grid(row=0, column=1, sticky="e")

        # --- Canvas area ---
        canvas_frame = tk.Frame(master, bg="#1e1e1e")
        canvas_frame.grid(row=1, column=0, sticky="nsew", padx=12, pady=6)

        self.canvas = tk.Canvas(canvas_frame, bg="#1e1e1e", highlightthickness=0, bd=0)
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<Configure>", self._on_canvas_resize)

        # --- Class buttons: three semantic rows ---
        self.btn_frame = tk.Frame(master, bg="#1e1e1e")
        self.btn_frame.grid(row=2, column=0, sticky="ew", padx=12, pady=(6, 6))

        row1 = ["E", "E/S0", "S0", "S0/a", "S"]                # basic
        row2 = ["SB0", "SB0/a", "SB", "SBa", "SBb", "SBc"]     # bars
        row3 = ["Sa", "Sb", "Sc", "Irr/Pec", "Star/other"]     # detailed/other

        maxcols = max(len(row1), len(row2), len(row3))

        lbl_classes = tk.Label(self.btn_frame, text="Select a galaxy class:",
                               font=self.font, fg="#cfcfcf", bg="#1e1e1e", anchor="w")
        lbl_classes.grid(row=0, column=0, columnspan=maxcols, sticky="w", pady=(0, 4))

        for c in range(maxcols):
            self.btn_frame.grid_columnconfigure(c, weight=1, uniform="cls")

        self.buttons = []

        def add_row(items, r):
            for c, label in enumerate(items):
                btn = tk.Button(self.btn_frame, text=label, font=self.font,
                                bg="#333333", fg="#ffffff",
                                activebackground="#555555", activeforeground="#ffffff",
                                relief="flat", bd=2, highlightthickness=0,
                                command=lambda lbl=label: self.classify(lbl))
                btn.grid(row=r, column=c, padx=4, pady=4, sticky="nsew")
                btn.bind("<Enter>", self.on_enter)
                btn.bind("<Leave>", self.on_leave)
                self.buttons.append(btn)

        add_row(row1, 1)
        add_row(row2, 2)
        add_row(row3, 3)

        for r in (1, 2, 3):
            self.btn_frame.grid_rowconfigure(r, pad=2)

        # --- Control buttons (separate row; no overlap) ---
        self.ctrl_frame = tk.Frame(master, bg="#1e1e1e")
        self.ctrl_frame.grid(row=3, column=0, sticky="ew", padx=12, pady=(4, 12))
        for c in range(3):
            self.ctrl_frame.grid_columnconfigure(c, weight=1, uniform="ctrls")

        self.Button_Undo = tk.Button(self.ctrl_frame, text="Undo", font=self.font,
                                     bg="#444444", fg="#ffffff",
                                     activebackground="#666666", activeforeground="#ffffff",
                                     relief="flat", bd=2, highlightthickness=0,
                                     command=self.undo_last)
        self.Button_Skip = tk.Button(self.ctrl_frame, text="Skip", font=self.font,
                                     bg="#444444", fg="#ffffff",
                                     activebackground="#666666", activeforeground="#ffffff",
                                     relief="flat", bd=2, highlightthickness=0,
                                     command=self.skip_image)
        self.Button_Quit = tk.Button(self.ctrl_frame, text="Exit", font=self.font,
                                     bg="#444444", fg="#ffffff",
                                     activebackground="#666666", activeforeground="#ffffff",
                                     relief="flat", bd=2, highlightthickness=0,
                                     command=self.ask_quit)

        self.Button_Undo.grid(row=0, column=0, padx=4, pady=4, sticky="nsew")
        self.Button_Skip.grid(row=0, column=1, padx=4, pady=4, sticky="nsew")
        self.Button_Quit.grid(row=0, column=2, padx=4, pady=4, sticky="nsew")

        for b in (self.Button_Undo, self.Button_Skip, self.Button_Quit):
            b.bind("<Enter>", self.on_enter)
            b.bind("<Leave>", self.on_leave)

        # Shortcuts
        self.master.bind("<BackSpace>", lambda e: self.undo_last())
        self.master.bind("<space>", lambda e: self.skip_image())
        self.master.bind("<Escape>", lambda e: self.ask_quit())

        # First image
        if self.index < self.n:
            self.update_image_display()
        else:
            tkMessageBox.showinfo("Info", "")
            self.master.destroy()

    # ---------- Logic ----------

    def update_image_display(self):
        path = self.all_paths[self.index]
        filename = os.path.basename(path)
        self.ImgTitle.config(text=filename)
        self.lbl_progress.config(text=f"{self.index + 1} / {self.n}")
        self.display_image(path)

    def display_image(self, path):
        try:
            pil = Image.open(path).convert("RGB")
            self._current_img_pil = pil
            self._draw_image_to_canvas(pil)
        except Exception as e:
            tkMessageBox.showerror("Error", f"Failed to load image:\n{e}")
            self.skip_image()

    def _on_canvas_resize(self, event):
        if self._current_img_pil is not None:
            self._draw_image_to_canvas(self._current_img_pil)

    def _draw_image_to_canvas(self, pil_img):
        # Canvas size
        cw = max(1, self.canvas.winfo_width())
        ch = max(1, self.canvas.winfo_height())

        # Uniform scale (never upscale)
        sx = cw / pil_img.width
        sy = ch / pil_img.height
        scale = min(1.0, sx, sy)

        img = pil_img
        if scale < 1.0:
            new_w = max(1, int(pil_img.width * scale))
            new_h = max(1, int(pil_img.height * scale))
            img = pil_img.resize((new_w, new_h), Image.Resampling.LANCZOS)

        # Center
        x0 = (cw - img.width) // 2
        y0 = (ch - img.height) // 2

        self._photo = ImageTk.PhotoImage(img)
        self.canvas.delete("all")
        self.canvas.create_image(x0, y0, anchor="nw", image=self._photo)

    def classify(self, selected_class):
        if self.index >= self.n:
            self._finalize()
            return
        name = self.all_paths[self.index]
        self.map_cls[name] = selected_class
        self._write_tsv()
        self.index += 1
        if self.index < self.n:
            self.update_image_display()
        else:
            self._finalize()

    def skip_image(self):
        if self.index < self.n:
            self.index += 1
            if self.index < self.n:
                self.update_image_display()
            else:
                self._finalize()

    def undo_last(self):
        if self.index > 0:
            self.index -= 1
            name = self.all_paths[self.index]
            if name in self.map_cls:
                del self.map_cls[name]
                self._write_tsv()
            self.update_image_display()

    def _write_tsv(self):
        # Write only entries in order of all_paths
        tmpfile = self.pngfile_out.with_suffix(self.pngfile_out.suffix + ".tmp")
        rows = [(p, self.map_cls[p]) for p in self.all_paths if p in self.map_cls]
        df = pd.DataFrame(rows, columns=['path',' class'])
        df.to_csv(tmpfile, sep=',', index=False, header=self.include_header)
        os.replace(tmpfile, self.pngfile_out)

    def _finalize(self):
        self._write_tsv()
        tkMessageBox.showinfo('Classification completed', f'Data saved to: {self.pngfile_out}')
        self.master.destroy()

    def ask_quit(self):
        if tkMessageBox.askokcancel("Exit", "Are you sure you want to exit?\nData will be saved."):
            self._write_tsv()
            self.master.destroy()

    def on_enter(self, event): event.widget.configure(bg="#666666")
    def on_leave(self, event): event.widget.configure(bg="#333333")


if __name__ == "__main__":
    # Example:
    # cgPNG("./imagenes_sdss", "class2.txt", include_header=False)
    pass
