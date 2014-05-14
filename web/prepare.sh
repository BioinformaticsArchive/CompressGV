#!/bin/bash

sed -i "s|GIT-COMMIT-HASH|`git log --pretty=format:'%h' -n 1`|" index.html;
analytics=`cat analytics.js`;
analytics=`printf '%q' $analytics`;
sed -i "s|GOOGLE-ANALYTICS|${analytics}|" index.html;

jsMD5=`md5sum CompressGV.js | awk '{print $1}'`;
find -maxdepth 1 -iname 'CompressGV_md5_*.js' | xargs -r rm;
cp CompressGV.js "CompressGV_md5_${jsMD5}.js";
sed -i "s|CompressGV.js|CompressGV_md5_${jsMD5}.js|" index.html;

cp index.html index_raw.html;
indexMD5=`md5sum index_raw.html | awk '{print $1}'`;
sed -i "s|data-version='master'|data-version='${indexMD5}'|" index.html;
