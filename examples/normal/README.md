zmat2xyz ethanol.zmat
zmat2xyz spce.zmat 
packmol < pack.inp 
create_conf -m 2.0 simbox.xyz 360 ethanol.xyz oplsaa.ff 360 spce.xyz spce.ff
mpirun -np 4 /opt/lammps/lmp_openmpi -i in.lmp
