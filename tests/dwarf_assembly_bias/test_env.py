import os, dwarf_assembly_bias

def test_env():
    print('Python PATH:', os.environ['PYTHONPATH'])    
    
def test_version():
    print(dwarf_assembly_bias.version)