#From http://www.blossomassociates.com/molbio/

#awk -F "\t" -f ~/lib/revcomp.awk --source '{print revcomp($0)}'

function complement(b) {
    if ( "c" == b ) b = "g";
    else if ( "g" == b ) b = "c";
    else if ( "a" == b ) b = "t";
    else if ( "t" == b ) b = "a";
    else if ( "u" == b ) b = "a";

    else if ( "C" == b ) b = "G";
    else if ( "G" == b ) b = "C";
    else if ( "A" == b ) b = "T";
    else if ( "T" == b ) b = "A";
    else if ( "U" == b ) b = "A";

    else if ( "m" == b ) b = "k";
    else if ( "r" == b ) b = "y";
#   else if ( "w" == b ) b = "w";
#   else if ( "s" == b ) b = "s";
    else if ( "y" == b ) b = "r";
    else if ( "k" == b ) b = "m";
    else if ( "v" == b ) b = "b";
    else if ( "h" == b ) b = "d";
    else if ( "d" == b ) b = "h";
    else if ( "b" == b ) b = "v";
#   else if ( "x" == b ) b = "x";
    else if ( "n" == b ) b = "x";

    else if ( "M" == b ) b = "K";
    else if ( "R" == b ) b = "Y";
#   else if ( "W" == b ) b = "W";
#   else if ( "S" == b ) b = "S";
    else if ( "Y" == b ) b = "R";
    else if ( "K" == b ) b = "M";
    else if ( "V" == b ) b = "B";
    else if ( "H" == b ) b = "D";
    else if ( "D" == b ) b = "H";
    else if ( "B" == b ) b = "V";
#   else if ( "X" == b ) b = "X";
    else if ( "N" == b ) b = "X";

#   else if ( "." == b ) b = ".";

    return(b);
}

# Return the reverse complement of a sequence of bases.
function revcomp( theBases ) {
  answer = "";
  l = length( theBases );
  for ( i = l; 0 < i; i-- ) {
    b = substr( theBases, i, 1 );
    b = complement(b);

    answer = answer b;
  }
  return answer;
}


# Return the reverse of a sequence of bases.
function rev( theBases ) {
  answer = "";
  l = length( theBases );
  for ( i = l; 0 < i; i-- ) {
    b = substr( theBases, i, 1 );

    answer = answer b;
  }
  return answer;
}
