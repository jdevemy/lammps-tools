    zmat2xyz CC.zmat
    zmat2xyz CO.zmat
    zmat2xyz spce.zmat
    fusion_mols CC.xyz CO.xyz CC-CO.ff
    packmol < pack.inp
    create_conf --freetype fep01 -m 2.0 simbox.xyz 1 CC.xyz/CO.xyz CC-CO.ff 360 spce.xyz spce.ff
    mpirun -np 4 /opt/lammps/lmp_openmpi -i in.lmp
