;;; basilisk-mode.el

;;; Copyright: (C) 2018 Stephane Popinet
;; 
;;     This program is free software; you can redistribute it and/or
;;     modify it under the terms of the GNU General Public License as
;;     published by the Free Software Foundation; either version 2 of
;;     the License, or (at your option) any later version.
;;     
;;     This program is distributed in the hope that it will be useful,
;;     but WITHOUT ANY WARRANTY; without even the implied warranty of
;;     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
;;     GNU General Public License for more details.
;;     
;;     You should have received a copy of the GNU General Public License
;;     along with GNU Emacs; if not, write to the Free Software
;;     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
;;     02110-1301 USA
;;
;; To use this package, you can save this file somewhere in your
;; load-path and put the following in your .emacs at a minimum:
;;
;; (require 'basilisk-mode)

(require 'markdown-mode)
(require 'mmm-mode)

(define-derived-mode basilisk-mode c-mode "Basilisk C"
  "Major mode for editing Basilisk C files."
  (mmm-mode)
)

;; customize mmm-mode
(mmm-add-classes
 '((basilisk
    :submode markdown-mode
    :face mmm-declaration-submode-face
    :front "^[ \t]*/[*][*]"
    :back "[*]/[ \t]*$"
      :include-front t
      :include-back t
      )))
(mmm-add-mode-ext-class 'basilisk-mode 'nil 'basilisk)

;; customize markdown-mode
(setq markdown-command "page2html")
(setq markdown-command-needs-filename t)
(add-hook 'markdown-mode-hook
	  `(lambda ()
	     (local-set-key [f8]  'markdown-export-and-preview)
	     (local-set-key [f9]  'markdown-export)
	     ))

(provide 'basilisk-mode)
