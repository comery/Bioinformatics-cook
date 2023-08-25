#!/usr/bin/perl

#创建1
%fruit = ("apples",17,"bananas",9,"oranges","none");

#创建2
%fruit = ("apples"=>17,"bananas"=>9,"oranges"=>"none");

#创建3
@fruit = ("apples",17,"bananas",9,"oranges","none");
%fruit = @fruit;


$fruit{"lime"} = 1; #增加
delete$fruit{"lime"}; #删除

@fruitsubs = keys(%fruits); #取出所有键
@fruitindexes = sort keys(%fruits); #取出所有键，并按字母排序
@fruitvalues = values(%fruits); #取出所有值


# 方式一
foreach $holder (keys(%records)){
	$record = $records{$holder};
	#操作
} 

#方式二
while (($holder, $record) = each(%records)) {
	#操作
}


#创建
@array = (0) x 5;
@array = 5..9;
@array = ('a', 'b');
@array = qw/9 0 2 1 0/;


push(@array, $i); #从数组的末尾加入元素$i
$i = pop(@array); #从数组的末尾取出元素，存入$i
$i = shift(@array);  #从数组的开头取出元素，存入$i
unshift(@array, $i); #从数组的开头加入元素$i

@array = split(/\t/, $str); #将字符串$str以\t为分割符切分，存入数组@array
$str = join("\t", @array);  #将数组@array以\t为连接符连接，存入字符串$str

scalar(@array); # 总数目
$#array;        # 最大下标，总数目-1

@arrayo = ("let", "less", "much");
@changes = map { uc($_) } @arrayo; #遍历并执行{}中的操作，这里是变大写uc
@selects = grep { /l/ } @arrayo;   #遍历并按{}中操作筛选，这里是匹配字母l
print "@changes\n";
print "@selects\n";


@array = ( 0 .. 6 );
@array1 = ( 'a' .. 'd' );
@replaced = splice( @array, 3, 2, @array1 );
print "replaced:     @replaced\n",
      "with:         @array1\n",
      "resulting in: @array\n\n";




@array = (8, 2, 32, 1, 4, 16);
print join(' ', sort { $a <=> $b } @array), "\n"; #升序
print join(' ', sort { $b <=> $a } @array), "\n"; #降序

sub numerically { $a <=> $b };
print join(' ', sort numerically @array), "\n";  #函数

@languages = qw(fortran lisp c c++ Perl python java);
print join(' ', sort { $a cmp $b } @languages), "\n"; # ASCII升序
print join(' ', sort { $b cmp $a } @languages), "\n"; # ASCII降序

$str = join "\n", 
	map { $_->[0] }
		sort { $a->[1] <=> $b->[1] }
			map { [ $_, (split)[-1] ] }
				split /\n/, $str; 

