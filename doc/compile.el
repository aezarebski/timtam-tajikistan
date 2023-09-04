(defun compile-and-clean (base-filename)
  "Compile LaTeX document with pdflatex and bibtex,
   clean auxiliary files, and make a backup of the output."

  ;; Compile LaTeX document
  (shell-command (format "pdflatex %s" base-filename))
  (shell-command (format "bibtex %s" base-filename))
  (shell-command (format "pdflatex %s" base-filename))
  (shell-command (format "pdflatex %s" base-filename))

  ;; List of extensions of temporary files to remove
  (setq aux-extensions '("aux" "log" "out" "bbl" "blg" "run.xml"))

  ;; Remove temporary files
  (dolist (ext aux-extensions)
    (let ((tmp-file (format "%s.%s" base-filename ext)))
      (when (file-exists-p tmp-file)
        (delete-file tmp-file))))

  ;; Make a backup copy of the pdf with today's date
  (copy-file (format "%s.pdf" base-filename)
             (format "%s-%s.pdf" base-filename (format-time-string "%Y-%m-%d"))))

(defun export-and-compile (org-file-path)
  "Export an Org-mode file to LaTeX and then compile it."
  (with-current-buffer (find-file-noselect org-file-path)
    ;; Export to LaTeX
    (org-latex-export-to-latex)
    (let ((base-filename (file-name-sans-extension org-file-path)))
      ;; Call compile-and-clean on the LaTeX file
      (compile-and-clean base-filename))))

(export-and-compile "~/projects/tajikistan-timtam/doc/readme.org")
