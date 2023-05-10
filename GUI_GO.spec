# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(
    ['GUI_GO.py'],
    pathex=[],
    binaries=[],
    datas=[('.\\global_optimization.py', '.'), ('.\\data_structure.py','.'),('.\\material_structure.py','.'),
    ('.\\material_model.py','.'), ('.\\default_script.txt','.'), ('.\\form_factor.pkl','.'),
    ('.\\form_factor_magnetic.pkl','.'), ('.\\Perovskite_Density.txt','.'), ('.\\Atomic_Mass.txt','.'),
    ('.\\demo.h5','.'), ('.\\User_Guide_v0.3.pdf','.'), ('.\\license.txt','.')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='GUI_GO',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='GUI_GO',
)
