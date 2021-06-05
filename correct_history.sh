#!/bin/sh
git filter-branch --env-filter --force '
OLD_EMAIL1="vinicius@biotonix.com"
OLD_EMAIL2="36831295+vbordalo@users.noreply.github.com"
CORRECT_NAME="vbordalo"
CORRECT_EMAIL="vinicius.bordalo@gmal.com"
if [ "$GIT_COMMITTER_EMAIL" = "$OLD_EMAIL1" ] or [ "$GIT_COMMITTER_EMAIL" = "$OLD_EMAIL2" ]
then
    export GIT_COMMITTER_NAME="$CORRECT_NAME"
    export GIT_COMMITTER_EMAIL="$CORRECT_EMAIL"
fi
if [ "$GIT_AUTHOR_EMAIL" = "$OLD_EMAIL" ] or [ "$GIT_AUTHOR_EMAIL" = "$OLD_EMAIL2" ]
then
    export GIT_AUTHOR_NAME="$CORRECT_NAME"
    export GIT_AUTHOR_EMAIL="$CORRECT_EMAIL"
fi
' --tag-name-filter cat -- --branches --tags
