
set grid
set title 'Bin central, fichero 20 ns'
set xlabel 't (s)'

dt=20e-6 #s
T1=1000e-3 #s
T2=200e-3 #s
t10=2000*dt
A1=1-0.04484
t20=2000*dt
A2=0.87729

f1(t)=(1-A1*exp(-(t-t10)/T1))
f2(t)=A2*exp(-(t-t20)/T2)
f(t)=sqrt(f1(t)*f1(t) + f2(t)*f2(t))

p 'output_ix1iy1t_gvec.txt' u ($1*dt):(sqrt($3*$3+$4*$4+$5*$5)) w lp t '|v|'
rep 'output_ix1iy1t_gvec.txt' u ($1*dt):3 w l t 'vx'
rep 'output_ix1iy1t_gvec.txt' u ($1*dt):4 w l t 'vy'
rep 'output_ix1iy1t_gvec.txt' u ($1*dt):5 w lp t 'vz'
rep 'output_ix1iy1t_gvec.txt' u ($1*dt):(sqrt($3*$3+$4*$4)) w lp t '|vxy|'
rep f1(x) lw 3 t 'f1(t)=exp(-t/T1)'
rep f2(x) lw 3 lc 12 t 'f2(t)=1-exp(-t/T2)'
rep f(x) lw 3 t 'sqrt(f1**2 + f2**2)'
