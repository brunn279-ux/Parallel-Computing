if has("gui_running")
    set guioptions-=T    " Turn off toolbar in gvim
 endif

set nocompatible     " Only use vim not vi
"set hidden           " Hide buffers instead of closing them
"set number           " Turn on line numbers
set ls=2             " Always show status line
set ruler            " Show the cursor position all the time
set visualbell t_vb= " Turn off error beep/flash
set title            " Show title in console title bar
set modeline         " Last lines in document sets vim mode
set modelines=5      " Number lines checked for modelines
"set shortmess=atI    " Abbreviate messages
"set wrap             " Wrap text lines
"set linebreak        " Break at specific char instead of width char
set nolist           " List disables linebreak
set tabstop=8        " Number of spacer per tab (affects existing tabs)
set shiftwidth=3     " Tab spaces on reindent
set softtabstop=8    " How many columns for new tabs
set expandtab        " Replace tabs with spaces
"set history=1000     " Number of commands remembered
set undolevels=1000  " yep
"set exrc             " Allow remote vimrc
set secure           " Only execute autocommand calls if owned by user
"set autoindent       " Always set autoindenting on
"set smartindent      " Smart indent
set backspace=indent,eol,start
set textwidth=100

" Some useful vim alterations
 let mapleader=","    " use the comma for executing plugin commands
 syntax on            " Turn on syntax highlighting
" Use ; instead of shift ; = :, e.g. :w (4 keystrokes), ;w  (2 keystrokes)
nnoremap ; :
" select text in between delimiters when using %
"noremap % v%

"     If you want to do more extensive vim customization,
"     install the vundle vim pacakge manager like this:
"       mkdir -p ~/.vim/bundle
"       cd ~/.vim/bundle
"       git clone https://github.com/VundleVim/Vundle.vim.git
"       cd
"     You can then uncomment the following lines, selecting the
"     Plugin commands for the packages you want to try, or adding
"     other packages that interest you.  When you have uncommented
"     or added what you like, run this from the command line to
"     download and install all the plugins you've requested:
"       vim +PluginInstall +qall

" Don't highlight the words that you have searched for.
set nohlsearch

" A function for automatically inserting include guards into new header files.
 function! Insert_gates()
    let gatename = substitute(toupper(expand("%:t")), "\\.", "_", "g")
    execute "normal ggi#ifndef " . gatename . "__\n#define " . gatename . "__\n"
    execute "normal Go\n#endif // " . gatename . "__"
    normal kk
 endfunction


set formatoptions=ql

" Expand tabs into spaces for most files.
set expandtab

" Enable file type detection.
" Use the default filetype settings, so that mail gets 'tw' set to 72,
" 'cindent' is on in C files, etc.
" Also load indent files, to automatically do language-dependent indenting.
filetype plugin indent on

" Filetypes
" (au=autocmd,setl=setlocal,ts=tabstop,sw=shiftwidth,sts=softabstop)
au FileType helpfile set nonumber " no line numbers when viewing help

autocmd FileType markdown,text,c,cpp,rst setlocal spell spelllang=en_us
autocmd FileType c,cpp  set formatoptions=ql cindent comments=sr:/*,mb:*,el:*/,://
autocmd BufNewFile *.{h,hpp,hh} call Insert_gates()

au BufNewFile,BufReadPost *.cpp,*.hpp,.h,*.cc,*.hh set filetype=cpp
au BufNewFile,BufReadPost *.md set filetype=markdown

au FileType python setl expandtab ts=4 sw=4 sts=4
au FileType python nnoremap <buffer> <F7> :exec '!python' shellescape(@%, 1)<cr>
au FileType python noremap <silent> <F8> :!clear;python %<CR>
"au FileType python setl omnifunc=syntaxcomplete#Complete

au Filetype gitcommit setl spell textwidth=72
autocmd BufNewFile,BufRead *.cls  setf tex
autocmd FileType cls  setf tex
" Turn on spelling for LaTeX
autocmd FileType tex setlocal spell spelllang=en_us

" Allow tabs for temporary files (like the ones created by CVS and Perforce).
autocmd BufNewFile,BufRead /tmp/*,/private/tmp/*,/var/tmp/* set noexpandtab

"let g:markdown_fenced_languages = ['html', 'python', 'bash=sh', 'c', 'cpp']

set list
set encoding=utf-8
set listchars=tab:‣◦,trail:◇
set pastetoggle=<F2>
" Change indent options to make constructors better
set cinoptions=(0,u0
" Change indent options to make constructors better
 let c_comment_strings=0

set bg=dark
" Maybe these are for light backgrounds?
"hi DiffChange ctermbg=DarkMagenta guibg=DarkMagenta ctermfg=white guifg=white
"hi DiffText ctermbg=Red guibg=Red ctermfg=white guifg=white
"hi DiffAdd ctermbg=DarkBlue guibg=DarkBlue ctermfg=white guifg=white
hi clear SpellCap
hi clear SpellBad
"hi SpellCap term=reverse ctermfg=red ctermbg=white  gui=undercurl guisp=NONE
"hi SpellBad term=reverse ctermfg=red ctermbg=yellow  gui=undercurl guisp=NONE
hi SpellCap term=underline cterm=underline gui=undercurl guisp=NONE
hi SpellBad term=underline cterm=underline gui=undercurl guisp=NONE
hi clear SpellLocal

autocmd BufNewFile,BufRead *.git/COMMIT_EDITMSG    setf gitcommit
autocmd BufNewFile,BufRead *.git/config,.gitconfig setf gitconfig
autocmd BufNewFile,BufRead git-rebase-todo         setf gitrebase
autocmd BufNewFile,BufRead .msg.[0-9]*
     \ if getline(1) =~ '^From.*# This line is ignored.$' |
     \   setf gitsendemail |
     \ endif
autocmd BufNewFile,BufRead *.git/**
     \ if getline(1) =~ '^\x\{40\}\>\|^ref: ' |
     \   setf git |
     \ endif

set wildignore=*.d,*.o

set laststatus=2

map \gq ?^$\\|^\s*\(\\begin\\|\\end\\|\\item\\|\\label\)?1<CR>gq//-1<CR>
omap lp ?^$\\|^\s*\(\\begin\\|\\end\\|\\item\\|\\label\)?1<CR>//-1<CR>.<CR>

map <F4> :%s/\s\+$//e<CR>
