This assumes that you have successfully generated initial conditions and you are setting up an entirely new set up initial conditions that haven't been ran before.  IF you are only adding any more bins to already existing initial conditions, please see the end of this document!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

What one needs to run, you need to wait for one command to finish (they run quickly) before running the next:


./setup.sh  a1  a2  a3 a4 a5
./setup2.sh a1 a2 a3 a4 a5 b6 b7 b8

cd v-USPhydro2

./writeGEN.sh c1 a2 c3 0 999

That's all. The last command will be reran multiple times enough you've ran over all bins.  That one you can get multiple bins up and running at the same time.

Now for what that means...

a1 Ion SNN directory folder name
a2 Varied parameters e.g. EOS21
a3- START subfolders
a4- END subfolders
a5- Grid step size

b6- eta/s
b7- fitting constant
b8- output options 0 (normal), 3 print off all energy profiles, 1=QM fluctuations

c1- ions name up to number i.e. XeXe=XeXe5440TeVdef_b236
HOWEVER, if you are running different deformations, rename the ions! NOTE, they must be all with a-z c1 can't handle #'s i.e.
UtUt=UtUt200GeV, then c1=UtUt and a1=UtUt200GeV.  Otherwise, you can't run different deformations at the same time!

c3- Bin number to run (you can only set up so many events to run at a time)

An Example (you will need to change names according to the initial conditions):

./setup.sh XeXe5440TeVdef_b236 EOS21 0 10 0.05

./setupIC2.sh XeXe5440TeVdef_b236 EOS21 0 10 0.05 0.047 119 0

cd v-USPhydro2

./writeGEN.sh XeXedef EOS21 0 0 999
./writeGEN.sh XeXedef EOS21 1 0 999
...
./writeGEN.sh XeXedef EOS21 10 0 999


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JUST ADDING BINS TO ALREADY EXISTING INITIAL CONDITIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

./moveIC.sh a1 a2 a3 a4 a5

cd v-USPhydro2

./writeGEN.sh c1 a2 c3 0 999
