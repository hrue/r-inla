#!/usr/bin/perl

## copied from git@github.com:clearlinux/make-fmv-patch.git

# Copyright (c) 2016, Intel Corporation

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# This program will make its best effort on finding the candidates sources for FMV
# patch them and create a bunch of ready to apply patches.

# USAGE: perl make-fmv-patch.pl <path_to_logfile> <full_path_to_source_code>

my %f;
my $fmv = {};
my ($log_file, $source_path) = @ARGV;

sub patch_function {

    my $attribute     = '__attribute__((target_clones(INLA_CLONE_TARGETS "default")))' . "\n";
    my $head          = '#define INLA_CLONE_TARGETS "avx2","arch=atom",' . "\n";
    my $target_clones = 'targets.conf';

    if (-e $target_clones) {
        open(my $fh, '<:encoding(UTF-8)', $target_clones) or die "Could not open file '$target_clones' $!";
        while (my $row = <$fh>) {
            chomp $row;
            $head = $head . $row . "\n";
        }
    }

    my ($file, @patch_line) = @_;
    my $patch_file = "$file" . "~";

    print "patching $file @ lines \(@patch_line\)\n";
    open(my $in, "<", $file) or die "$! - $file\n";
    open(my $out, ">", "$patch_file") or die __LINE__, " - $!\n";

    print $head;
    foreach (@patch_line) {
        my $line_num = $_;

        while (<$in>) {
            print $out $_;
            last if $. == ($line_num - 1);
        }

        print $out $attribute;
    }
    while (<$in>) {
        print $out $_;
    }

    close $out;
    close $in;

    my $diff = `diff -su $file $patch_file`;
    open(my $d, ">", "$file.patch") or die __LINE__, " - $!\n";
    print $d $diff;
    close($d);
    `rm $patch_file`;
}

sub find_file {
    my ($file) = @_;

    if ($file =~ /.*(\.\.\/)/p) {
	$file = ${^POSTMATCH};
    }

    my $n = (split('/',$file))[-1];
    chomp (my @matches = `find $source_path -iname $n`);

    foreach (@matches) {
        if ($_ =~ $n) {
            $file = $_;
        }
    }
    return $file;
}

open(BUILD_LOG, '<', "$log_file") or die $!;
while (<BUILD_LOG>) {
    if($_ =~ /(\S+):([0-9]*):([0-9]*): (optimized|note): (basic block|loop) (vectorized)/) {
	$fmv->{s_name} = (split('/',$1))[-1];
	$fmv->{f_name} = $1;
	push @{$f{$1}->{v_line}},$2;
    }
}
close(BUILD_LOG);

foreach (keys %f) {
    # print "$_ => @{$f{$_}->{v_line}}\n";
    my $fname = $_;
    my @flines = @{$f{$_}->{v_line}};

    if ($fname =~ /(\.c$)/ || $fname =~ /(\.cpp$)/) {
        my $f_path = &find_file($fname);

        my @keys = 0;
        my $i = 0;
        foreach (`ctags --c-kinds=f -x $f_path`) {
            chomp $_;
            $keys[$i++] = (split(/\s+/,$_))[2]; # get function number line
        }
        @keys = sort {$a <=> $b} @keys;

        my @match_line = 0;
        # foreach line vectorized look for its closest function
        foreach my $i ( 0 .. $#flines ) {
            foreach (@keys) {
            if ($_ < $flines[$i] ) {
                $match_line[$i] = $_;
            }
            }
        }

        my %h;
        foreach (@match_line) {$h{$_} = $_;} @match_line = keys %h;
        @match_line = sort {$a <=> $b} @match_line;

        &patch_function($f_path,@match_line);
    }
}
