#!/bin/sh

RM=$(git remote update 2>&1)
UPSTREAM=${1:-'@{u}'}
LOCAL=$(git rev-parse @)
REMOTE=$(git rev-parse "$UPSTREAM")
BASE=$(git merge-base @ "$UPSTREAM")

if [ $LOCAL = $REMOTE ]; then
    echo "Up-to-date"
elif [ $LOCAL = $BASE ]; then
    echo "Need to pull, please run git pull and recompile!"
    SHORT=$(git rev-parse  --short=10 HEAD)
    echo "Change log:"
    git log $SHORT\..origin/HEAD --oneline --decorate --color --abbrev=10
elif [ $REMOTE = $BASE ]; then
    echo "Need to push"
    SHORT=$(git rev-parse  --short=10 HEAD)
    echo "Change log:"
    git log $SHORT\..origin/HEAD --oneline --decorate --color --abbrev=10
else
    echo "Diverged"
fi
