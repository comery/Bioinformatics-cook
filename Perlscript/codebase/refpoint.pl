#!perl -w

my $ref_to_AoA = [
    [ "fred", "barney", "pebbles", "bamm bamm", "dino"],
    [ "homer", "bart", "marge", "maggie"],
    [ "george", "jane", "elroy", "judy"],
];

my @Aoa = (
      ["fred", "barney"],
      ["george", "jane", "elroy"],
      ["homer", "marge", "bart"],
      {"test"=>"aaa","test2"=>"bbb"}
);

print $ref_to_AoA->[2][3];  # 等价于 $ref_to_AoA->[2]->[3];
print $Aoa[3]{"test"};      # 等价于 $Aoa[3]->{"test"} 


