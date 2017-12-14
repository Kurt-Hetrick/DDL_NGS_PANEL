# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

umask 0007

export PATH=.:/bin:/usr/bin:$PATH

# User specific aliases and functions
module load gcc/5.1.0
module load sge/2011.11p1

#######################################
##### .bash_history file settings #####
#######################################

# Make the .bash_history file very big
export HISTSIZE=1000000
export HISTFILESIZE=1000000

# appends session history to .bash_history
shopt -s histappend

# immediately add commands to .bash_history instead of on logout.
# This is helpful since I multiple concurrent sessions running.
PROMPT_COMMAND='history -a'

# I've seen the extended version of above...but am not implementing it at this time.
# I think this will re-read .bash_history to all open terminals...
# After each command, append to the history file and reread it
# export PROMPT_COMMAND="${PROMPT_COMMAND:+$PROMPT_COMMAND$'\n'}history -a; history -c; history -r"

# Avoid duplicates in .bash_history
export HISTCONTROL=ignoredups:erasedups

# This sets the prompt
export PS1="\u@\h> "
