### External

The *external* directory contains external resources (libraries, programs, utilities, etc...) that are **not** part of the FOODIE library, but that are used into the tests suite facility. Each external resource has its own license, thus FOODIE users are kindly encouraged to read them.

#### External resources

The current external resources list contains:

+ [FLAP](https://github.com/szaghi/FLAP):  *Fortran command Line Arguments Parser for poor people*; it is integrated into the FOODIE project as a *git submodule* (it being completely independent from FOODIE) thus it can be easily update within FOODIE, see the following paragraph;
+ [IR_Precision](https://github.com/szaghi/IR_Precision): *Pure Fortran (2003+) library for ensuring codes portability*; it is integrated into the FOODIE project as a *git submodule* (it being completely independent from FOODIE) thus it can be easily update within FOODIE, see the following paragraph;
+ [pyplot-fortran](https://github.com/jacobwilliams/pyplot-fortran): *For generating plots from Fortran using Python's matplotlib.pyplot*; it is integrated into the FOODIE project as a *git submodule* (it being completely independent from FOODIE) thus it can be easily update within FOODIE, see the following paragraph;
+ [WenOOF](https://github.com/Fortran-FOSS-Programmers/WenOOF): *WENO interpolation Object Oriented Fortran library*; it is integrated into the FOODIE project as ai regular *git clone* of its mainstream (it being **not** completely independent from FOODIE) thus it cannot be easily update within FOODIE: to update it a manual *git pull* is necessary, see the following paragraph;

#### Update the external resources

For the *git submodule* external resources the update is done with
```shell
git submodule update --init --recursive
```
executed into the project root.

For the *git cloned* external resources the update is done with
```shell
cd foodie-root/external/cloned-resource/
git pull
```
where *foodie-root* is the local root of FOODIE project and *cloned-resource* is the actual resource root.
