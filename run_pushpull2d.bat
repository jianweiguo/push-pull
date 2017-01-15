set PNUMBER=1024
set DMIN=0.85
set RC=0.67
set SDA=0.02
set ITER=3500

FOR /L %%s in (0,1,0) DO (
pushpull2d.exe -p tmp/ -I 0 -i %ITER% -d %DMIN% -r %RC% -q 012 -R %%s -met %PNUMBER%
)

pause