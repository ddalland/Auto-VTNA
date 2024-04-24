# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['Auto_VTNA_GUI.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=['openpyxl.cell._writer', 'scipy.optimize.curve_fit', 'scipy.sparse'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['torch', 'keras', 'tensorflow', 'NumExpr', 'torchvision', 'conda', 'setuptools', 'psutil', 'pytest', 'pygments', 'pyarrow', 'IPython', 'h5py', 'bokeh', 'sqlalchemy', 'lz4', 'tables', 'PyYAML', 'traitlets', 'anyio', 'lxml', 'botocore', 'dask', 'qtpy', 'cryptography', 'bcrypt', 'argon2', 'jinja2', 'nacl', 'notebook', 'boto3', 'numpydoc', 'sphinx', 'babel', 'sqlite3', 'nbformat', 'jsonschema', 'win32com', 'pythoncom'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='Auto_VTNA_GUI',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
