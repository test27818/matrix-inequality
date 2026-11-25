"""Check original finite-kernel endpoint then complete signed-curvature replay.
Use: python3 reconstructed/run_signed_curvature.py EDGE FIRST_CAP BPU Q_LEVELS LEVELS
C++20 is needed for this faster implementation; binary is kept in a temporary directory.
"""
import os,subprocess,sys,tempfile
from pathlib import Path
from verify_signed_curvature import endpoint,P
HERE=Path(__file__).resolve().parent
def main():
    edge,cap,bpu,nq,levels=map(int,sys.argv[1:]) if len(sys.argv)>1 else (4,4,8,16,6)
    assert 1<=cap<=edge<=8 and 1<=bpu<=128 and (bpu&(bpu-1))==0 and 1<=nq<=1024 and (nq&(nq-1))==0 and 1<=levels<=10
    atom=endpoint(edge,1,max(32,bpu))
    num=atom*(1<<40)//P
    print('scaled cutoff-atom lower='+str(num)+'/1099511627776',flush=True)
    with tempfile.TemporaryDirectory() as folder:
        binary=str(Path(folder)/'signed-replay')
        subprocess.run([os.environ.get('CXX','c++'),'-std=c++20','-O3','-fsanitize=undefined','-fno-sanitize-recover=all',str(HERE/'signed_curvature_fast.cpp'),'-o',binary],check=True)
        return subprocess.run([binary,*map(str,(edge,cap,bpu,nq,levels,num))]).returncode
if __name__=='__main__':raise SystemExit(main())
