set history=1000
set savehist=100


# -------------------------------------------------------
# for AFNI: auto-inserted by init_user_dotfiles.py

# add AFNI abin to PATH
setenv PATH ${PATH}:$HOME/abin

# set up tab completion for AFNI programs
# (only do this in an interactive shell)
if ( $?prompt ) then
   if ( "$prompt" != "" ) then
      if ( -f $HOME/.afni/help/all_progs.COMP ) then
         source $HOME/.afni/help/all_progs.COMP
      endif
   endif
endif
# -------------------------------------------------------

setenv R_LIBS $HOME/sw/R-4.3.1
