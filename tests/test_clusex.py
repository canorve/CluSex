

import pytest
import subprocess as sp


from clusex import pysex 




# simple run test
def test_exit():
    with pytest.raises(SystemExit) as e:
        pysex.main()
    assert e.type == SystemExit 
    assert e.value.code == 2 


def test_sextractor():

    runcmd="sex"
    err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  

    assert err.returncode == 0, "is Sextractor installed?"


def test_ds9():

    runcmd="ds9 -quit -lower" 
    err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  
    assert err.returncode == 0, "is Ds9 installed?"




