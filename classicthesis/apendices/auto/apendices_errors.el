;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "apendices_errors"
 (lambda ()
   (LaTeX-add-labels
    "appex_errors"
    "app_errors:eq:taui"
    "eq:WF"
    "eq:tail"
    "eq:Stau"))
 :latex)

