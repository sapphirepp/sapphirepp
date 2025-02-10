# IDE {#ide}

@tableofcontents

In our experience the complexity of developing @sapphire is strongly reduced if
an integrated development environment (IDE) is used. An IDE can be thought of as
a software suite that comprises tools to handle the different tasks appearing
when writing software. In particular, an IDE allows users to easily navigate the
many source files a project may exist of, integrates version control systems
like [git](https://git-scm.com/), unifies the format of the code and provides
users with a front-end to a debugger, e.g. [gdb](https://sourceware.org/gdb/).
As the navigation, the auto-completion and the formatting relies on a language
server, e.g. [clangd](https://clangd.llvm.org/), a modern IDE needs to be able
to communicate with it using the [language server
protocol](https://en.wikipedia.org/wiki/Language_Server_Protocol) (lsp).

## Emacs

There are various IDEs available that fulfil these requirements. A very popular
choice is [VS Code](https://code.visualstudio.com/). [Eclipse
IDE](https://eclipseide.org/) is also widely used. I chose GNU [Editor
Macros](https://www.gnu.org/software/emacs/) (Emacs). Its design differs much
from VS Code and Eclispe, see the Section [General
Architecture](https://en.wikipedia.org/wiki/Emacs#General_architecture) of the
Wikipedia page for a concise discussion. Prejudices about Emacs being
out-of-date or dead, are prejudices. Emacs does everything other IDEs do. The
difference is that Emacs follows a Do-It-Yourself (DIY) philosophy and in this
way it provides users with an extreme level of flexibility. A user has the
ability to fine tune his programming tool to perfectly suite his or her needs.
Moreover, using Emacs means the keyboard is key. I find not having to use the
mouse very convenient. Of course, GNU Emacs is free and open-source software and
I experience the community around it as very vivid, friendly and inspiring.

That said, the Emacs setup I present here is always preliminary and shows the
essential parts of my personal setup. It is actually a minimal setup that shall
provide the reader with a starting point from which she may explore the plethora
of options Emacs and its packages offer. Moreover, Emacs is developing fast and
I change my Emacs setup constantly to keep up with the new reflect this
development. I will update this page whenever important of my Emacs
configuration change. A reader familiar with Emacs will notice, that the
presented setup is conservative in the sense that I mainly use packages which
are ported with the current Emacs version. However, I try to mention alternative
choices whenever I know some. If you think, that I miss an essential package or
if you have a better setup, please write me an
[email](mailto:nils.schween@mpi-hd.mpg.de) and share it with me.

### Basics

For the presented configuration to work, at least Emacs 29.1. is needed. Emacs
is configured with a `.emacs` or a `init.el` file and the presented
[elisp](https://en.wikipedia.org/wiki/Emacs_Lisp) snippets need to copied to one
of them.

The basic configuration concerns all modes and variables that can be set without
loading an additional packages.

```elisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Emacs configuration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package emacs
  :config
  ;; add the melpa archive to the list of package archives
  (add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/"))

  ;; automatically inserts the second parenthesis
  (electric-pair-mode)

  ;; highlights the matching parenthesis when point is before opening or after closing one
  (show-paren-mode)

  ;; wrap at 80 characters
  (setq-default fill-column 80)

  ;; sources of the treesitter language grammars
  (setq treesit-language-source-alist
	'((cpp "https://github.com/tree-sitter/tree-sitter-cpp")
	  (c "https://github.com/tree-sitter/tree-sitter-c")
	  (cmake "https://github.com/uyha/tree-sitter-cmake")
	  (python "https://github.com/tree-sitter/tree-sitter-python")
	  ))

  ;; Remap major modes to their corresponding treesitter-mode
  (setq major-mode-remap-alist
	'((c-mode . c-ts-mode)
	  (c++-mode . c++-ts-mode)
	  (cmake-mode . cmake-ts-mode)
	  ))
  )
```

Package configuration is done itself with the package `use-package`. Type `C-h R
use-package` to browse its manual. Packages are either part of Emacs (in this
case they are listed in [ELPA](https://elpa.gnu.org/)), or they can be
downloaded from package archives like the [NonGNU
ELPA](https://elpa.nongnu.org/) or the [MELPA](https://melpa.org/#/). You can
get a list of all available packages with `M-x list-packages`.

The last two blocks are needed to use
[Tree-sitter](https://tree-sitter.github.io/tree-sitter/) for the parsing of
code. Emacs traditionally relied for this on regular expressions. Tree-sitter
unifies the interface for all programming modes and simplifies their
implementation.

```elisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; prog-mode configuration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package prog-mode
  :hook ((prog-mode . display-fill-column-indicator-mode)
	 (prog-mode . display-line-numbers-mode)
	 (prog-mode . hl-line-mode)
	 (prog-mode . highlight-indentation-current-column-mode)
	 ))

```

The major programming modes, e.g. `python-mode`, are very often derived from `prog-mode`. This means
that all settings done in the `prog-mode` block apply to them.

### Completion

Completion may heavily speed-up typing. Many packages implement completion UIs
for the mini-buffer. Particular popular choices are
[helm](https://emacs-helm.github.io/helm/) and
[https://oremacs.com/swiper/](ivy). A minimalistic option is
[vertico](https://github.com/minad/vertico), which I am currently using. I would
like to mention that I used ivy + swiper before and that I found it to be very
reliable as well. If vertico is used, it is worth to take a look at the
[consult](https://github.com/minad/consult) package. It provides search
and navigation commands that feed vertico.

In-buffer completion can, for example, be achieved with the
[company](https://company-mode.github.io/) or the
[corfu](https://github.com/minad/corfu) package.

All packages are highly customizable and it is worth to read their documentation
if a specific behaviour is (un-)desired.

```elisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Vertico - VERTical Interactive COmpletion
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package vertico
  :init
  (vertico-mode)
  :ensure t

  ;; Different scroll margin
  ;; (setq vertico-scroll-margin 0)

  ;; Show more candidates
  ;; (setq vertico-count 20)

  ;; Grow and shrink the Vertico minibuffer
  ;; (setq vertico-resize t)

  ;; Optionally enable cycling for `vertico-next' and `vertico-previous'.
  ;; (setq vertico-cycle t)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Corfu - COmpletion in Region FUnction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package corfu
  ;; Optional customizations
  ;; :custom
  ;; (corfu-cycle t)                ;; Enable cycling for `corfu-next/previous'
  ;; (corfu-quit-at-boundary nil)   ;; Never quit at completion boundary
  ;; (corfu-quit-no-match nil)      ;; Never quit, even if there is no match
  ;; (corfu-preview-current nil)    ;; Disable current candidate preview
  ;; (corfu-preselect 'prompt)      ;; Preselect the prompt
  ;; (corfu-on-exact-match nil)     ;; Configure handling of exact matches

  ;; Enable Corfu only for certain modes. See also `global-corfu-modes'.
  ;; :hook ((prog-mode . corfu-mode)
  ;;        (shell-mode . corfu-mode)
  ;;        (eshell-mode . corfu-mode))

  ;; Recommended: Enable Corfu globally.  This is recommended since Dabbrev can
  ;; be used globally (M-/).  See also the customization variable
  ;; `global-corfu-modes' to exclude certain modes.
  :init
  (global-corfu-mode))
```


<div class="section_buttons">

| Previous                                  |                                  Next |
|:------------------------------------------|--------------------------------------:|
| [Coding Convections](#coding-conventions) | [Acknowledgements](#acknowledgements) |

</div>

@author Nils Schween (<nils.schween@mpi-hd.mpg.de>)
@date Mo 10. Feb 11:53:46 CET 2025
