# SelfForce-1D

SelfForce-1D is a code infrastructure for simulating Extreme Mass Ratio Inspirals using the effective source approach to the self-force problem. Currently, only a scalar charge in a Schwarzschild spacetime background has been implemented, but it is the hope that more systems will be added soon.

This is a complete rewrite of 1DScalarWaveDG. The new code is much more modular in order to make it easy to extend the code with other systems of equations.

These set of files also contain new python files like logpage_2, store_2 ad parser_2. These files are modified versions of logpage, store and parser found in the Einstein Toolkit. These files make it possible to create HTML pages that contain test results. The directory structure (records, docs, etc) is also the same as the Einstein Toolkit.

## SelfForce-1D documentation

For more info, visit https://bitbucket.org/peterdiener/selfforce-1d/src/master/.

## GitHub Actions & files for CI testing

  - [Introduction](#introduction)
    - [What is Github Actions](#what-is-github-actions)
    - [Why Github Actions](#why-github-actions)
    - [File Overview](#file-overview)
  - [Explanation of Files](#explanation-of-files)
    - [parser_2](#parser_2)
    - [store_2](#store_2)
    - [logpage_2](#logpage_2)
  
## Introduction

### What is Github Actions

Github Actions is continuous integration/ continous development platform that runs
a set of commands on a repository. Github Actions allows the creation of user 
created modules that automates certain commonly used workflows. On each push, the 
workflow is run in a docker container (running ubuntu).

### Why Github Actions

- Github Actions allows tests to be run on their servers as such there is no server maintenance required
- There is less security risks because its hosted on the cloud rather than an active server
- Easier local testing allowing for new features to be tested easier
- Flexibility to to tailor the reports to the Einstein Toolkit since we can design our own parsers and tools.
- Larger community giving more opportunity for more plugins than Jenkins.


### File Overview

- `parser_2.py` - Parses the log files
- `logpage_@.py` - Generates the HTML pages
- `store_2.py` - Stores logs for future use
- `mail.py` - Send email each time tests are run
- `new_test_nums.csv` - Stores summary stats from logs (note:currently has random input fed to generate sample of bokeh plot)
- `records/` - Folder contains compilation logs, logs with summary of tests, and individual test logs and diffs. 
- `docs/new_index.html` - HTML page that is displayed on mojamil.github.io/einsteintoolkit/

## Explanation of Files

### parser_2

This python script is used to parse the log files for required data.
The parsing is done by using regex expression matching to find the
necessary information. A brief description of what each function
does

`create_summary(file)` This function looks for the summary of the tests stored in log files such
as build__2_1.log or build__1_2.log:

![summary](https://github.com/Hrishikesh-Kalyanaraman/einsteintoolkit/blob/SelfForce-1D/images/log_file_sample.png)

`get_tests(file)` Gets the name of the test that passed and failed as listed in log files such
as build__2_1.log:

![pass-fail](https://github.com/Hrishikesh-Kalyanaraman/einsteintoolkit/blob/SelfForce-1D/images/pass-fail.png)

`test_comp(file1,file2)` Compares the passed and failed tests from file1 and file2 and returns
which files are newly passing,newly failing, newly added and removed.

`get_data(file)` Retrieves singular field of data from a csv and returns it as a list

### store_2

`copy_tests(test_dir,version,procs)`  function that is currently commented. Need link to logs anf diffs files of selfforce-1d to continue (can currently implement by linking to .out files)

`copy_logs(version)` This copies the test logs for future use

`copy_compile_log(version)` This copies the compilation logs for future use

`copy_index(version)`  This copies the old html files showing test results for future use

`get_version()` Gets the version based on the stored files if there are no stored files
returns 1

### logpage_2

Logpage.py generates tables for the html report page and outputs as an html page as
shown here:

https://rhaas80.github.io/einsteintoolkit/

This file gets the last few commits using githubs REST API for commits and workflow runs as 
shown in these documentation links: https://docs.github.com/en/rest/reference/repos#commits and https://docs.github.com/en/rest/reference/actions

This file uses bokeh, a python library, to generate plots. The plots are created using python code and bokeh
then converts to javascript and html.

![bokeh](https://github.com/mojamil/einsteintoolkit/blob/gh-pages/images/bokeh.png)
![plot](https://github.com/mojamil/einsteintoolkit/blob/gh-pages/images/plot.PNG)

Bokeh's plotting works similar to other plotting libraries. First a figure is generated and attributes can
be added such as tools to zoom, labels, axis ranges, etc. Bokeh plots using glyphs i.e. given data it will
plot it in the format specified for example p.line shown above generates a line graph and p.circle can be
used for scatter plots. Bokeh can show its plot locally and save it as a file or generate html and javascript
for the plot as shown below:

![bokeh2](https://github.com/mojamil/einsteintoolkit/blob/gh-pages/images/bokeh2.png)
