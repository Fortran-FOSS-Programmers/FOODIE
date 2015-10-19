### How-to perform automated errors analysis

In order to automate the errors analysis computation a bash script is provided, the `errors_analysis.sh` script. The scripts automates the following steps:

1. build the tests programs and place the executable into the errors analysis directory;
1. run the tests with the auto-errors-analysis feature enabled (currently for only the oscillation test);
1. save a summary log of the analysis;
1. if Tecplot program is present into PATH, the results of the errors analysis are plotted by means of the provided layouts and exported into PNG images in order to be easily used into the LaTeX manuscript under writing.

To run the script from the FOODIE root directly do:

```shell
cd papers/announcement/errors_analysis/
./errors_analysis.sh
```

For exporting PNG images by means of Tecplot some auxiliary scripts and macro are provided. They are placed into `papers/announcement/errors_analysis/utilities`. The results of the errors analysis are placed into the newly created directories named accordingly to the test program name, namely:

1. papers/announcement/errors_analysis/results-burgers
1. papers/announcement/errors_analysis/results-euler-1D
1. papers/announcement/errors_analysis/results-lorenz
1. papers/announcement/errors_analysis/results-oscillation

The PNG images into the results directories are sym-linked into the `papers/announcement/images/errors-analysis/` sub directories.
