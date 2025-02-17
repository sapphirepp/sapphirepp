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
from VS Code and Eclipse, see the Section [General
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

### Basics {#basics}

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

  ;; COrfu
  ;; Enable indentation+completion using the TAB key.
  ;; `completion-at-point' is often bound to M-TAB.
  (setq tab-always-indent 'complete)

  ;; Tree-sitter
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

### Completion

Minibuffer completion may heavily speed-up choosing canidates when opening files
or executing commands. Many packages implement completion UIs for the
minibuffer. Particular popular choices are
[helm](https://emacs-helm.github.io/helm/) and
[https://oremacs.com/swiper/](ivy). A minimalist option is
[vertico](https://github.com/minad/vertico), which I am currently using. I would
like to mention that I used ivy + swiper before and that I found it to be very
reliable as well. If vertico is used, it is worth to take a look at vertico's
extensions and the [consult](https://github.com/minad/consult) package. It
provides search and navigation commands that feed vertico.

In-buffer completion can, for example, be achieved with the
[company](https://company-mode.github.io/) or the
[corfu](https://github.com/minad/corfu) package. It helps to keep track of
complicated variable names.

I note that all mentioned packages are highly customisable and it is worth to
read their documentation if a specific behaviour is (un-)desired.

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

To explore the new completion UI you can, for example, type `M-x` and execute an
arbitrary command. In the [Basics](#basics) section, we changed the value of the
variable `tab-always-indent` to `complete`. This allows to a spawn a child frame
with possible completion candidates using `TAB`. If you prefer auto-completion,
please read the
[Auto-completion](https://elpa.gnu.org/packages/doc/corfu.html#Auto-completion)
section in corfu's manual.

### Programming

The major programming modes, e.g. `python-mode`, are very often derived from `prog-mode`. This means
that all settings done in the `prog-mode` block apply to them.

```elisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; prog-mode configuration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package prog-mode
  :hook (
     ;; Activates the fill-column indication display, i.e. a vertical line is
	 ;; drawn at the end of the column set by `fill-column` 
	 (prog-mode . display-fill-column-indicator-mode)
	 ;; Displays the line numbers at the left border of the screen
	 (prog-mode . display-line-numbers-mode)
	 ;; Highlights the current line 
	 (prog-mode . hl-line-mode)
	 ;; Adds visible columns to ease the recognition of the indentation level
	 (prog-mode . highlight-indentation-current-column-mode)
	 ;; Activate spell checking of comments in code
	 (prog-mode . flyspell-prog-mode)
	 ))

```
The spell checking of the comments while coding is done with the package
[flyspell](https://www.gnu.org/software/emacs/manual/html_node/emacs/Spelling.html).
I note that flyspell only works if a spell checker, e.g.
[Aspell](http://aspell.net/), and a corresponding dictionary are installed. If
your operating system has the [Enchant](https://rrthomas.github.io/enchant/)
library installed, you can consider the spell checker
[Nuspell](https://nuspell.github.io/) in conjunction with the Emacs package
[jinx](https://github.com/minad/jinx). 

Nested parentheses are very common. The package
[rainbow-delimiters](https://github.com/Fanael/rainbow-delimiters) helps to know
how far I went down in the hierarchy of parentheses. I note that
rainbow-delimiters is a MELPA package.

```
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Rainbow Delimiters -  have delimiters be colored by their depth
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package rainbow-delimiters
  :ensure t
  :hook (prog-mode . rainbow-delimiters-mode))
```

Eglot is the ported Emacs client for the “Language Server Protocol” (LSP). The name
“Eglot” is an acronym that stands for “Emacs polyGLOT”.

```
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Eglot: Interface to language servers LSP
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package eglot
  :bind (("C-x l" . eglot)
     :map eglot-mode-map
     ("C-x l w r" . eglot-reconnect)
     ("C-x l w q" . eglot-shutdown)
     ("C-x l w k" . eglot-shutdown-all)
     ("C-x l = r" . eglot-format)
     ("C-x l = =" . eglot-format-buffer)
     ("C-x l r r" . eglot-rename)
     ("C-x l a a" . eglot-code-actions)))
```

A widely used and more feature rich alternative to eglot is
[lsp-mode](https://emacs-lsp.github.io/lsp-mode/). The cumbersome shortcuts for
the eglot commands are chosen to be close to the keymap used by lsp-mode. I used
lsp-mode for quite some time and I was very satisfied. However, since Emacs 29.1
eglot is part of the ported emacs packages and I decided to switch, since I only use a
small subset of options language servers offer and eglot seems to cover them.
For a LSP client to work, a language server is necessary to which the client can
connect.

For C/C++ I am using the [clangd](https://clangd.llvm.org/) server. An
alternative option is [ccls](https://github.com/MaskRay/ccls). If installed,
clangd should be automatically detected by eglot. If you have ccls and clangd
installed, you can specify in your emacs configuration which one you would like to use with
the variable `eglot-server-programs`. If you like to start the language server
automatically whenever you open a C/C++ file, you need to add a hook to the corresponding
major modes, i.e. `c-ts-mode` and `c++-ts-mode`; see the section [Editor
plugins](https://clangd.llvm.org/installation#editor-plugins) of the clangd
website. Clangd needs to know the build flags of your software projects, see
[Project setup](https://clangd.llvm.org/installation#project-setup). If you use
CMake, as we do for @sapphire, you can create a `compile_commands.json` file,
that needs to be in the build directory of your project. 

For Python ...

To use eglot, open a source file in one of your (version controlled) software
projects and type `C-x l`. The source code will be indexed and eglot is ready to
use. Eglot functionality builds on other Emacs packages: 

- It shows the documentation of symbols (functions, variables etc.) at point using
[ElDoc](https://www.gnu.org/software/emacs/manual/html_node/emacs/Programming-Language-Doc.html),
- It detects coding errors and suggests fixes with the help of
[flymake.](https://www.gnu.org/software/emacs/manual/html_node/emacs/Flymake.html).
Type `C-x l a a` to apply the fixes.
- It allows users to jump to definitions of symbols at point with
[xref](https://www.gnu.org/software/emacs/manual/html_node/emacs/Xref.html) with
`M-,`. To jump back type `M-.`. It is not anymore necessary to generate TAGS.
- It allows buffer navigation by function, class, method etc. via
  [Imenu.](https://www.gnu.org/software/emacs/manual/html_node/emacs/Imenu.html).
  Type `M-g i` to jump to a function.
- It feeds corfu with completion candidates. 
- It formats your buffer with clang-format. Type `C-x l = =`.
- It allows you to easily rename symbol in your whole project, type `C-x l r r`

And more, see [Eglot features](https://www.gnu.org/software/emacs/manual/html_node/eglot/Eglot-Features.html).

Very often a software project does not exist of a single file. It has root
directory out of which a tree of folder files grow. In this case the
pre-installed
[project.el](https://www.gnu.org/software/emacs/manual/html_node/emacs/Projects.html)
comes very handy. To get you started, I list the commands I regularly use:

- `C-x p f` Visit a file that belongs to the current project
- `C-x p c` Run compilation in the current projects root directory
- `C-x p &` Run shell command asynchronously in the current project's root
directory 

For a complete list of commands, see [Project
Commands](https://www.gnu.org/software/emacs/manual/html_node/emacs/Project-File-Commands.html)
in the Emacs manual. I have to admit that I am missing a command to start a
debugger in the project's root directory. For an alternative project managment
tool, I refer the reader to the package [projectile](https://github.com/bbatsov/projectile). 

To get an overview of the files in a directory, you can use the already
installed package
[dired](https://www.gnu.org/software/emacs/manual/html_node/emacs/Dired.html).
If you prefer a tree view of your project, I suggest to use the package
[treemacs](https://github.com/Alexander-Miller/treemacs). 

Many people like to fold their code. Emacs comes the
[HideShow](https://www.gnu.org/software/emacs/manual/html_node/emacs/Hideshow.html)
mode to do this. If you prefer a code folding based using treesitter, take a
look at [treesit-fold](https://elpa.nongnu.org/nongnu/treesit-fold.html).

To debug my code I used [gdb](https://sourceware.org/gdb/). Emacs has a great
graphical interface to gdb. Just type `M-x gdb` and then `M-x gdb-many-windows`.
Breakpoints can be set by clicking into the fringe of your code windows. 

Peope working on clusters, sometimes like to use their emacs setup to remotely
edit code. This is possible with the package
[tramp](https://www.emacswiki.org/emacs/TrampMode) .

Very often, we have to write boiler plate code. The package
[https://github.com/joaotavora/yasnippet](yasnippet) offers a good amount of
code snippets. Though, I have to admit that I am not using them very often.

```
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; yasnippet
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package yasnippet
  :ensure t
  :init
  (yas-global-mode 1)
  )

(use-package yasnippet-snippets
	     :ensure t
	     :requires yasnippet
	     )
```

### Version control 

We use git to control the versions of our code. [Magit](https://magit.vc/) is
wonderful git interface. It does everything I need with only a few key strokes.
The package [forge](https://magit.vc/manual/forge.html) allows me to interact
with the API of [github](https://github.com), i.e. can open issues,
pull-requests etc. It also supports other git forges like
[gitlab](https://about.gitlab.com/). Read its manual to set it up.

```
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Magit - A git porcelain inside emacs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package magit
  :ensure t
  :bind ("C-x g" . magit-status)
  :config
  (setq magit-diff-refine-hunk t)
  )

(use-package forge
  :after magit
  :ensure t
  )
```

Open a file in a git controlled folder and type `C-x g` to open magit. If you
are familiar with git, magit is self-explanatory. Type `?` to open a transient
buffer that shows the git commands that you can call from magit. For example, a
simple work flow may be to edit a file, then stage the changes and subsequently
commit. In magit this boils to the key sequence `s` (staging), `c c` (commiting)
and subsequently entering a commit message. When you are done type `C-c C-c`. If
you like to your changes to a remote type `P p`.

With forge you can update your list of pull-requests, issues etc. by typing
sequence `N f f`.

### Documentation

We use markdown to write our documentation. Emacs has a very good markdown mode.
```
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;            Markdown-Mode          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package markdown-mode
  :ensure t
  :commands (markdown-mode gfm-mode)
  :mode (("README\\.md\\'" . gfm-mode)
         ("\\.md\\'" . markdown-mode)
         ("\\.markdown\\'" . markdown-mode))
  :init
  (setq markdown-command "multimarkdown")
  ;; Hashes only on the left-hand side of the heading
  (setq markdown-asymmetric-header t)
  ;; Font lock for inline and display LaTeX math expressions
  (setq markdown-enable-math t)
  )
```

### Auxiliary packages

In this section, I include packages that do not directly contribute to the IDE
capabilities of Emacs, but that proved to be useful in my daily programmer life.

An example, is the package `which-key`. Whenever you begin a series of key
strokes, it opens a new window at the bottom of of your Emacs frame, showing
possible continuations. This highly enhances the visibility of possible commands
provided by the various packages. It is automatically activated, just type `C-x`
to see its effect.

```
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; which-key
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(use-package which-key
  :ensure t
  :config
  (which-key-mode))
```

- jinx for spell checking comments
- multiple cursors

<div class="section_buttons">

| Previous                                  |                                  Next |
|:------------------------------------------|--------------------------------------:|
| [Coding Convections](#coding-conventions) | [Acknowledgements](#acknowledgements) |

</div>

@author Nils Schween (<nils.schween@mpi-hd.mpg.de>)
@date Mo 10. Feb 11:53:46 CET 2025
