
import inspect
import warnings

warnings.filterwarnings("ignore")


# Copyright 2025 Université Paris Cité, France
# Author: [Guillaume Rousseau](https://www.linkedin.com/in/grouss/), Physics Department, Paris, France
# This is supplemental materials and replication package associated with the preprint available on 
# - arXiv (https://arxiv.org/abs/2501.10145)
# - SSRN  (http://ssrn.com/abstract=5191689
# Current version of python scripts and ressources are available on [github's author page](https://github.com/grouss/growing-network-study)
# This work is currently licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)

def DisplayCopyrightInfo(Verbose=True):
    """
    This function display copyright and license information.
        Parameters :
            Verbose : Default=True. 
                      Print info
        Returns :
            outputstring :string
    """
    outputstring=80*"-"+"\n"
    
    outputstring+="Copyright 2025 Université Paris Cité, France \n"
    outputstring+="Author: Guillaume Rousseau, Physics Department, Paris, France \n"
    outputstring+="\n"
    outputstring+="(https://www.linkedin.com/in/grouss/)\n"
    outputstring+="\n"
    outputstring+="This archive contains the supplemental materials and replication package associated with the preprint available on :\n"
    outputstring+="- arXiv (https://arxiv.org/abs/2501.10145)\n"
    outputstring+="- SSRN  (http://ssrn.com/abstract=5191689\n"
    outputstring+="\n"
    outputstring+="Current version of python scripts and associated ressources are available on author's github page\n"
    outputstring+="(https://github.com/grouss/growing-network-study)\n"
    outputstring+="\n"
    outputstring+="This work is currently licensed under CC BY-NC-SA 4.0\n"
    outputstring+="(https://creativecommons.org/licenses/by-nc-sa/4.0)\n"
    outputstring+=80*"-"+"\n"
    if Verbose: 
        print(outputstring)
    else:
        return outputstring



def DisplayAllDoc(locallib):
    """
    Display the name and documentation of all user-defined functions within a given module.
    This utility filters and lists all functions defined in the provided module, 
    and prints their associated docstring, if available.
        Parameters :
            locallib : module
                       The module object to inspect.
    """
    functions_list = [func for func in dir(locallib) 
                  if inspect.isfunction(getattr(locallib, func)) 
                  and inspect.getsourcefile(getattr(locallib, func)) == locallib.__file__]
    print("Imported functions defined in",locallib.__file__,":")
    print("-"*80)
    for func_name in functions_list:
        func = getattr(locallib, func_name)
        doc = inspect.getdoc(func)
        print(f"Function: {func_name}()")
        if doc:
            print(f"{doc}")
        else:
            print("No documentation available")
        print("-"*80)
    
