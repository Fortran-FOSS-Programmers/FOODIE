---
project: FOODIE
src_dir: ../src/lib
src_dir: ../src/tests
output_dir: html/publish/
project_github: https://github.com/Fortran-FOSS-Programmers/FOODIE
project_download: https://github.com/Fortran-FOSS-Programmers/FOODIE/releases/latest
summary: Fortran Object oriented Ordinary Differential Equations integration library
author: Fortran-FOSS-Programmers Group
github: http://fortran-foss-programmers.github.io/
website: http://fortran-foss-programmers.github.io/
email: stefano.zaghi@gmail.com
md_extensions: markdown.extensions.toc(anchorlink=True)
               markdown.extensions.smarty(smart_quotes=False)
               markdown.extensions.extra
               markdown_checklist.extension
docmark: <
display: public
         protected
         private
source: true
warn: true
graph: true
sort: alpha
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!../README.md!}
---
