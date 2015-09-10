### External

The *external* directory contains external resources (libraries, programs, utilities, etc...) that are **not** part of the FOODiE library, but that are used into the tests suite facility. Each external resource has its own license, thus FOODiE users are kindly encouraged to read them.

#### External resources

The current external resources list contains:

+ [FLAP](https://github.com/szaghi/FLAP):  *Fortran command Line Arguments Parser for poor people*; it is integrated into the FOODiE project as a *git submodule* (it being completely independent from FOODiE) thus it can be easily update within FOODiE, see the following paragraph;
+ [IR_Precision](https://github.com/szaghi/IR_Precision): **; it is integrated into the FOODiE project as a *git submodule* (it being completely independent from FOODiE) thus it can be easily update within FOODiE, see the following paragraph;
+ [pyplot-fortran](https://github.com/jacobwilliams/pyplot-fortran): **; it is integrated into the FOODiE project as a *git submodule* (it being completely independent from FOODiE) thus it can be easily update within FOODiE, see the following paragraph;
+ [WenOOF](https://github.com/Fortran-FOSS-Programmers/WenOOF): **; it is integrated into the FOODiE project as ai regular *git clone* of its mainstream (it being **not** completely independent from FOODiE) thus it cannot be easily update within FOODiE: to update it a manual *git pull* is necessary, see the following paragraph;

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
where *foodie-root* is the local root of FOODiE project and *cloned-resource* is the actual resource root.
