domx=$1    
domy=$2    
domz=$3    
sizex=1
sizey=2
sizez=4      
res=$4        
codeTimeLim=$5
cycleLim=$6
tcool=0.1
Tmid=0.364
Tatm=0.636
zq=2.636
delta=2.211
R_AU=115.0
Sig_FUV=0.1
ionfrac_FUV=1.e-5
Am0=1.0
beta=1.e5


cat >> athinput.txt << EOF
<comment>
problem = Stratified 3D MRI with zero-net vertical flux
author  = Stone, Hawley, Gammie, & Balbus
journal =
config  = --with-problem=strat --with-eos=isothermal --enable-shearing-box

<job>
problem_id      = StratCooling      # problem ID: basename of output filenames
maxout          = 4          # Output blocks number from 1 -> maxout
num_domains     = 1          # number of Domains in Mesh

<output1>
out_fmt = hst                # History data dump
dt      = 62.831853

<output2>
out_fmt = tab                # tab data dump
dt      = 62.831853
out = prim
num = 0

<output3>
out_fmt = rst                # tab data dump
dt      = 62.831853               # time increment between outputs
num = 0

<output4>
name = 1d                
dt      = 0.62831853
num = 0


<time>
cour_no         = 0.4        # The Courant, Friedrichs, & Lewy (CFL) Number
nlim            = $cycleLim    # cycle limit
tlim            = $codeTimeLim      # time limit 

<domain1>
level           = 0         # refinement level this Domain (root=0)
Nx1             = $((sizex*res*2))       # Number of zones in X-direction
x1min           = -$sizex.0      # minimum value of X
x1max           = $sizex.0       # maximum value of X
bc_ix1          = 4         # boundary condition flag for inner-I (X1)
bc_ox1          = 4         # boundary condition flag for outer-I (X1)
NGrid_x1        = $domx         # with MPI, number of Grids in X1 coordinate

Nx2             = $((sizey*res*2))       # Number of zones in Y-direction
x2min           = -$sizey.0     # minimum value of Y
x2max           = $sizey.0       # maximum value of Y
bc_ix2          = 4         # boundary condition flag for inner-J (X2)
bc_ox2          = 4         # boundary condition flag for outer-J (X2)
NGrid_x2        = $domy         # with MPI, number of Grids in X2 coordinate

Nx3             = $((sizez*res*2))       # Number of zones in X3-direction
x3min           = -$sizez.0      # minimum value of X3
x3max           = $sizez.0       # maximum value of X3
bc_ix3          = 2         # boundary condition flag for inner-K (X3)
bc_ox3          = 2         # boundary condition flag for outer-K (X3)
NGrid_x3        = $domz         # with MPI, number of Grids in X3 coordinate

<problem>
omega           = 1.0
gamma           = 1.6666666666666667    # gamma = C_p/C_v
pres           = 0.5
dfloor         = 1.e-5
pfloor         = 1.e-6
tfloor         = 0.1

ipert           = 1          # 1 for random d,P, 2 for uniform Vx
amp             = 0.025      # dP/P <= amp

beta            = $beta       # Plasma beta
ifield          = 8          # 1 for zero-net-flux Bz, 2 for constant Bz
nwz             = 1

CASE            = 3          # Diffusivity: uniform (1) or user defined (2)
R_AU            = $R_AU
Sig_FUV         = $Sig_FUV
ionfrac_FUV     = $ionfrac_FUV
Am0             = $Am0


icool          = 2          # 0 for none, 1 for flat, 2 for cos profile
tauCool        = $tcool        # needs to be defined for any cooling
T0             = 0.5        # define for flat cooling
Tmid           = $Tmid        # define for cos cooling
Tatm           = $Tatm        # define for cos cooling
zq             = $zq        # define for cos cooling
delta          = $delta        # define for cos cooling




<log>
child_out_level=-1

EOF
