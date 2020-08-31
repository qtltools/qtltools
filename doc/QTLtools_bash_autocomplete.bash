#!/usr/bin/env bash
_QTLtools() 
{
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    mode="${COMP_WORDS[1]}"
    opts=$( QTLtools --help | grep -A9999999999 'Available modes:' | tail -n+2 | awk '{print $1}' | tr '\n' ' ' )
    if [[ ${prev} == QTLtools ]] ; then
    #if [[ ${COMP_CWORD} == 1 ]] ; then
        COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
        return 0
    elif [[ $mode != "" ]] && [[ " ${opts[@]} " =~ " ${mode} " ]] && [[ ${cur} == -* ]]; then
        local copts
        copts=$( QTLtools $mode --help | grep '^  \-' | awk '{if($3 ~ /^\-/) {print $1"\n"$3} else {print $1}}' | tr '\n' ' ')
        COMPREPLY=( $(compgen -W "${copts}" -- ${cur}) )
        return 0
    #else
    #    local copts
    #    copts=$(ls -a .)
    #    COMPREPLY=( $(compgen -W "${copts}" -- ${cur}) )
    #    return 0       
    fi
}
complete -o default -F _QTLtools QTLtools