# JuliaCompare.jl

JuliaCompare.jl is built on top of PromulaDBA.jl to simplify reading variables
from ENERGY 2020 DBA files.

## Installation

JuliaCompare.jl can be added to a project by running 

```julia
using Pkg; 
Pkg.add(url="https://github.com/pnvolkmar/JuliaCompare.jl")
```
## Usage

Before reading any variables from DBAs into Julia, the package must parse the 
source code from a specific ENERGY 2020 model. This is done with the `list_vars`
function which takes three inputs :

* the CODE_FOLDER, the Engine subdirectory where source files are stored
* the DATA_FOLDER, where DBAs are stored (unzipped)
* db_files, a list of files to be be parsed which create the DBA files in the DATA_FOLDER
  * NOTE: Most ENERGY 2020 Models use the same set of db_files which are listed in a vector called db_files that can be imported from JuliaCompare.

This function outputs a DataFrame `vars` listing each variable along with its dimensions, associated DBA file, and keys for each dimension (and where those keys can be found)

The `vars` DataFrame can be stored with a DATA_FOLDER as a `Loc` object. 

```julia
import JuliaCompare as J
import JuliaCompare: db_files
CODE_FOLDER = raw"path\to\Engine"
DATA_FOLDER = raw"path\to\DBA\files\directory"
vars = J.list_vars(CODE_FOLDER, DATA_FOLDER, db_files)
Loc1 = J.Loc(vars, DATA_FOLDER)
```

Once the `vars` DataFrame has been parsed and a `Loc` created, variables can be 
retrieved easily

```julia
TotPol = J.var("TotPol", Loc1)
xTotPol = J.var("xTotPol", Loc1)
diff = J.diff_var(TotPol, xTotPol)
J.plot_diff(diff)
```
