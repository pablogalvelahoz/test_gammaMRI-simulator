

set grid
set title 'transfer 30x30 20 ns'
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

ndots=100
column=2

p 'output_xyt_gvec.txt' u ($1*dt):(sqrt($87*$87+$88*$88+$89*$89)) every 1 w lines t '|v|'
#rep 'output_xyt_gvec.txt' u ($1*dt):87 every 1 w lines t 'vx'
#rep 'output_xyt_gvec.txt' u ($1*dt):88 every 1 w lines t 'vy'
rep 'output_xyt_gvec.txt' u ($1*dt):89 every 1 w lines t 'vz'
rep 'output_xyt_gvec.txt' u ($1*dt):(sqrt($87*$87+$88*$88)) every 1 w lines t '|vxy|'
#rep f1(x) lw 3 t 'f1(t)=exp(-t/T1)'
#rep f2(x) lw 3 lc 12 t 'f2(t)=1-exp(-t/T2)'
#rep f(x) lw 3 t 'sqrt(f1**2 + f2**2)'
rep 'output_xyt_gvec.txt' u ($1*dt):($86/200.) every 1 w l t 'gradient (no units)'

