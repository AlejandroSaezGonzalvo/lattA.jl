;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "ClassicThesis"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("scrreprt" "twoside" "openright" "titlepage" "numbers=noenddot" "headinclude" "footinclude" "cleardoublepage=empty" "abstract=on" "BCOR=5mm" "paper=a4" "fontsize=11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("silence" "") ("subcaption" "") ("graphicx" "") ("longtable" "") ("yfonts" "") ("adjustbox" "") ("rotating" "")))
   (TeX-run-style-hooks
    "latex2e"
    "classicthesis-config"
    "FrontBackmatter/Titlepage"
    "FrontBackmatter/Titleback"
    "FrontBackmatter/Dedication_2"
    "FrontBackmatter/Dedication"
    "FrontBackmatter/Abstract"
    "FrontBackmatter/Publications"
    "FrontBackmatter/Acknowledgments"
    "FrontBackmatter/Contents"
    "introduction/introduction"
    "introduction/introduccion"
    "cap2/lattice"
    "cap3/observables"
    "cap4/cap4"
    "cap5/cap5_new"
    "cap6/cap6_new"
    "conclusions/conclusions"
    "conclusions/conclusiones"
    "apendices/apendices_conventions"
    "apendices/apendices_SU3"
    "apendices/apendices_simulations"
    "apendices/apendices_solvers"
    "apendices/apendices_errors"
    "apendices/apendices_b"
    "apendices/apendices_GEVP"
    "apendices/apendices_a"
    "apendices/apendices_obs"
    "apendices/apendices_c"
    "apendices/apendices_f"
    "apendices/apendices_qm"
    "FrontBackmatter/Bibliography"
    "silence"
    "scrreprt"
    "scrreprt10"
    "subcaption"
    "graphicx"
    "longtable"
    "yfonts"
    "adjustbox"
    "rotating")
   (TeX-add-symbols
    "tmp"
    "oddsidemargin"
    "evensidemargin")
   (LaTeX-add-labels
    "pt:intro"
    "pt:fund"
    "pt:results"
    "pt:concl")
   (LaTeX-add-bibliographies
    "Bibliography"
    "publications"))
 :latex)

