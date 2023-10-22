# PROTACable
![Project Logo](./assets/logo/project-logo.png)
![Project Scheme](./assets/logo/project-scheme.png)

PROTACable is an end-to-end in-silico design toolkit for novel PROTACs

## Table of Contents

- [Preliminary](#preliminary)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Stage-1](#stage-1)
- [Stage-2](#stage-2)
- [Stage-3](#stage-3)
- [Stage-4](#stage-4)
- [Stage-5](#stage-5)
- [Contributing](#contributing)
- [References](#references)
- [License](#license)

## Preliminary

## Requirements

Different stages of PROTACable pipeline have variable set of requirements.
### Conda environment:
* See [Installation](#installation) below.
### Linux OS:
* This toolkit collection had been tested on Centos 8 OS.
* SLURM scheduler is needed for Stage 2 (protein-protein docking) due to the way the code was implemented.
### Software:
* Stage 1: GNINA v1.0 (optional).
* Stage 2: ProPOSE v2022.
* Stage 5: Maestro Schrödinger (optional).
### Hardware:
* GPU (RAM >=12GB) (optional)


## Installation

This toolkit collection had been tested on Centos 8 OS.
To get started with PROTACable, follow these instructions:

```
git clone https://github.com/giaguaro/PROTACable.git 
cd PROTACable
export PROTACable=${pwd}
```
> **Note:** This is a standout point.

(Note: Remove the extra backticks around the code blocks above.)

## Usage

Here's a basic example of how to use this project:

\```
from project_title import AwesomeFeature

# Initialize
feature = AwesomeFeature()

# Use the feature
result = feature.do_something_cool()
print(result)
\```

(Note: Remove the extra backticks around the code blocks above.)

## Features

### Feature 1

![Feature 1 Screenshot](https://path-to-your-image/feature1.png)

A brief description of this feature.

### Feature 2

Sub-headings can be used to describe sub-features or components:

#### Component 1
Details about this component.

#### Component 2
Details about this component.

## References

1. [Reference Title](https://www.example.com) - A short description of this reference.
2. Author Name, "Paper Title", Journal Name, Year. [Link to paper](https://www.example-paper.com)

## License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.
