# vertMixing_workFlow

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr

# Notes:

These notes are largely from [Getting started with workflowr][].

[Getting started with workflowr]: https://cran.r-project.org/web/packages/workflowr/vignettes/wflow-01-getting-started.html

- The two required subdirectories are `analysis/` and `docs/`. These directories should never be removed from the workflowr project.
- `analysis/`: This directory contains all the source R Markdown files for implementing the data analyses for your project. It also contains a special R Markdown file, index.Rmd, that does not contain any R code, but will be used to generate index.html, the homepage for your website.
- `docs/`: This directory contains all the HTML files for your website. The HTML files are built from the R Markdown files in `analysis/`. Furthermore, any figures created by the R Markdown files are saved here.

The optional directories are `data/`, `code/`, and `output/`. These directories are suggestions for organizing your data analysis project, but can be removed if you do not find them useful.

- `data/`: This directory is for raw data files.

- `code/`: This directory is for code that might not be appropriate to include in R Markdown format (e.g. for pre-processing the data, or for long-running code).

- `output/`: This directory is for processed data files and other outputs generated from the code and data. For example, scripts in `code/` that pre-process raw data files from `data/` should save the processed data files in `output/`.