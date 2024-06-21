# MOVIT-Source-Code

### Overview
This repository provides a demonstration on how to use the MOVIT `R` source code.

# System Requirements

## Software Requirements

### OS Requirements

The source code has been tested on Microsoft's Windows 10 operating system and Linux (Ubuntu 18.04). The source code should be compatible with Windows, Mac, and Linux operating systems.

Before using the MOVIT source code, users should have `R` version 4.3.0 or higher, and two packages installed.

### Installation  

First, we need to install two dependencies `ivreg` and `mvnfast`:

    install.packages(c("ivreg", "mvnfast"))
    
which should install within a couple of minutes on a standard machine.

# Demo

We first load all the source code dependencies:

```
library(ivreg)
library(mvnfast)
library(psych)
library(remMap)
#library(mvnfast) #for replicating simulation II only
```

and the source code containing all the main functions:

```
source("MOVIT.R")
```
