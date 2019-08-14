#!/bin/sh
# -*- coding: utf-8 -*-

#create the directory with the date of the file and move the file to this directory

shopt -s nullglob  # this line is so that it does not compain when no logfiles are found
for filename in MOD021*.hdf; do # Files considered are the ones startign with test and ending in .log
    foldername=$(echo "$filename" | awk '{print (substr($0, 19, 4));}'); # The foldername is characters 19 to 23 from the filename (if they exist)
    mkdir -p "$foldername"  # -p so that we dont get "folder exists" warning
    mv "$filename" "$foldername"
    echo "$filename $foldername" ;
done

for filename in MOD02*.hdf; do # 
	for foldername in ????; do
		filepart=$(echo "$filename" | awk '{print (substr($0, 16, 4));}'); 
		if [ "$foldername" == "$filepart" ]; then
			mv "$filename" "$foldername"
			echo "$filename $foldername" ;
		fi
	done
done

for filename in MOD03*.hdf; do # 
	for foldername in ????; do
		filepart=$(echo "$filename" | awk '{print (substr($0, 16, 4));}'); 
		if [ "$foldername" == "$filepart" ]; then
			mv "$filename" "$foldername"
			echo "$filename $foldername" ;
		fi
	done
done
