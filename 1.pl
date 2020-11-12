#!/usr/bin/perl -w
# -*-mode:cperl -*-
package main;

=pod
Description: Code for https://www.fi.muni.cz/~hanka/ai/du_1.pdf
Usage: difficult.. everything hardcoded
Date created: 2020-11
Author: Vojtech Prusa
=cut
use strict;
use warnings;
# TODO:
use Array::Utils qw(:all);
use Data::Dumper qw(Dumper);

use 5.016003; # current version, TODO update ...
########################################################################
package SLog;
use Data::Dumper qw(Dumper);

use constant {
  ALG                  => 'Alg',
  ALG_EVO              => 'ALG_EVO',
  ALG_EVO_LOOP         => 'ALG_EVO_LOOP',
  ALG_EVO_LOOP_ISSIUF  => 'EVO_ISSUF',
  # ALG_EVO_LOOP_ISSUF   => 'ALG_EVO_LOOP_ISSUF',

  ALG_EVAL_GSWN        => 'ALG_EVAL_GSWN',
  ALG_EVAL_GSWN_VN     => 'ALG_EVAL_GSWN_VN',
  ALG_EVAL_GSWN_SUM    => 'ALG_EVAL_GSWN_SUM',
  ALG_EVAL_GSWN_SUM_VN => 'ALG_EVAL_GSWN_SUM_VN',
  # => '',
};

# TODO: Here comment anyone line you do not want to log
# Commenting 'L1b_L2A' but keeping 'L1b' will still log 'L1b_L2A' because of its prefix
# current setup should Dump:
=pod
TODO
=cut

my @SHOULD_LOG = (
  # ALG,
  ALG_EVO_LOOP,
  ALG_EVO_LOOP_ISSIUF,
  ALG_EVAL_GSWN,
  # ALG_EVAL_GSWN_VN,
  # ALG_EVAL_GSWN_SUM_VN,
  # ALG_EVAL_GSWN_SUM,
);

sub Log {
  if (scalar(@_) > 1) {
    my ($name, @val) = @_;
    my @prefix = split("_", $name);
    if ($name ~~ @SHOULD_LOG || $prefix[0] ~~ @SHOULD_LOG) {
      print "$name ";
      print Dumper(@val);
    }
    else {
    }
  }
  else {
    print Dumper(@_);
  }
}

package main;

# Wrapper for case when you do not want to call SLog::Log() but Log() because type speed matters
sub Log {
  SLog::Log(@_);
}

sub main {
  # https://www.fi.muni.cz/~hanka/ai/prusvitky/treti.pdf
  # ::7
  # TODO:
  # - its possible that i do smth bad with memory management, I was trying to keep as less memory as possible though...
  # --> review & performance analyses needed
  evalEvolutionAlg();
}

##################################################
# data structs,
# alg's variable's names has prefix 'd'
##################################################

# nodeX => {nodeY=> d(nodeX, nodeY), ...}, ...
# this results in duplicity of edge's weights
# Note!: Eh, index from 1...
# TODO ...
# another way to store data as triplets where first 2 are vertices and third is weight
# as in https://github.com/vprusa/MUNI-MA015-graph-algorithms/blob/e510fbe87b1c8d73f817fb13995bc3e3007f1ab3/algorithms/mst/directed/EdmondsBranching/EdmondsBranching.py#L13
# g.e. Edges: [(1,2,4),(1,3,5),(1,4,3),(2,3,3),(2,4,5),(3,4,4)]
# another another way is to use some strucutre to hold weight and just link it to given weight..
# current way is convenient for fast data lookup, although that could be done with previous methods + caching or smth...
my %dG = (
  1 => { 2 => 4, 3 => 5, 4 => 3 },
  2 => { 1 => 4, 3 => 3, 4 => 5 },
  3 => { 1 => 5, 2 => 3, 4 => 4 },
  4 => { 1 => 3, 2 => 5, 3 => 4 }
);

my @dWi = (10, 5, 1, 5);
my $dN = 4;
my $dV = 4;
my $dM = 1;


# TODO pass arguments from cmd line
my %settings = (
  isTest                    => 0, # test disables initial randomization

  reproductionCrossbreeding => 1,  # enables crossbreeding
  reproductionMutation      => 1,  # enables random mutations
  mutationRandMax           => 10, # TODO one var in % , see also mutationRandThreshold
  mutationRandThreshold     => 8,  # TODO one var in % , see also mutationRandMax

  maxIteration              => 100,
  maxSameResCnt             => 10
);

my %bestResSoFar = ();
my $iteration = 0;
my $sameResCnt = 0;

# choose initial solution
# my $rS=1; # for debug purposes, because using rand() for initial debug is not a good idea...
# my $rS=int(rand((keys %dG)-1))+1;
# Heuristics could rely on graph type (topology or other information)
# or for general graph i could consider some heuristic like
# - prioritizing/choosing nodes according to
# -- max(/min) number of neighbours
# -- lowest edge's weight
# -> Note: these heuristics that generate metadata from graph should be done in O(n) or max O(n^2), but not worse than
#  the rest of the algorithm itself...

##################################################
# subroutines
##################################################

my $getCostWithMultiNodesI = 0; # helper for debugging

sub getCostWithMultiNodes {
  Log(SLog::ALG_EVAL_GSWN, "getCostWithMultiNodes");
  ++$getCostWithMultiNodesI; # 13 -> res.totalCost=20

  # my ($limit, %genes) = @_;
  my ($limit, $genesRef) = @_;
  my %genes = %{$genesRef};

  my %res = ();
  $res{'genes'} = \%genes;

  my %visitedNodes = preparePreVisited();
  my $totalCost = 0;
  for (1 .. $dM) {
    my @gK = keys %dG;

    my @visitedNodesEmptyI = ();
    foreach my $n (keys %visitedNodes) {
      if ($visitedNodes{$n} != 0) {
        push(@visitedNodesEmptyI, $n);
      }
    }
    my @vnK = @visitedNodesEmptyI;
    my @unusedNodesA = array_diff(@gK, @vnK);
    my %unusedNodesH;
    $unusedNodesH{$_} = $dG{$_} for (@unusedNodesA); # just filter and cpy for now

    # my $nodeI = pickFreeNodeIndex(%unusedNodesH);
    # my $nodeI = pickFreeNodeIndex(%unusedNodesH);
    my $nodeI = 0;
    if ($_ == 1) {
      if ($genes{'startPos'} eq 'random') {
        $nodeI = randIndex(%unusedNodesH);
      }
      else {
        $nodeI = $genes{'startPos'};
      }
    }
    else {
      if ($genes{'nextPosOffset'} eq 'random') {
        $nodeI = randIndex(%unusedNodesH) % scalar @unusedNodesA;
      }
      else {
        $nodeI = $genes{'nextPosOffset'} % scalar @unusedNodesA;
      }
    }
    my ($newCost, %newVisitedNodes) = getCostWithNodes($nodeI, $limit, %visitedNodes);
    %visitedNodes = (%visitedNodes, %newVisitedNodes);
    if ($_ == 1) { # just once
      $newCost += $dWi[$nodeI - 1];
    }
    $totalCost += $newCost;
    # $res{$_} = ( 'cost' => $newCost, 'visNodes' => {%visitedNodes} );
    # $res{$_} = ( 'cost' => $newCost, 'visNodes' => {%visitedNodes} );
    $res{$_} = { 'cost' => $newCost, 'visNodes' => { %visitedNodes } };
  }
  $res{'totalCost'} = $totalCost;
  return %res;
}

sub getCostWithNodes {
  Log(SLog::ALG_EVAL_GSWN, "getCostWithNodes");
  my ($nodeI, $limit, %preVisited) = @_;

  my %node = %{$dG{$nodeI}};
  my $neighbourCost = 0;
  my %curVisited = (%preVisited);
  # TODO: should or should not I use 1 or 2 foreach loops sequentially?
  #  Because it may result in breaking on already visited nodes - "$preVisited{$neighbour} == 0" in case of loop in recursion
  #  Do I miss something? The simple example is insufficient for this Use Case, it would need some with possible loops..

  # idk if i need to copy the '%preVisited', maybe just '$curVisited{$nodeI}++;' now and '$curVisited{$nodeI}--;' later?
  $curVisited{$nodeI}++;
  # also add curent cost, TODO rename variable

  # TODO sort is quite expensive -> parametrization...
  foreach my $neighbourI (sort (keys %node)) {
    my $neighbourDistance = $node{$neighbourI};
    if ($preVisited{$neighbourI} == 0 && $limit >= $neighbourDistance) {
      $neighbourCost += $dWi[$neighbourI - 1];
      # %visitedNodes = (%visitedNodes, ($nodeI));
      $curVisited{$nodeI} = 1;
      my $newLimit = $limit - $node{$neighbourI};
      if ($newLimit > 0) {
        my ($newCost, %newVisitedNodes) = getCostWithNodes($neighbourI, $newLimit, %curVisited);
        $neighbourCost += $newCost;
        %curVisited = (%curVisited, %newVisitedNodes);
      }
      elsif ($newLimit == 0) {
        # TODO brainstorm if propagation of curVisited should not be somewhere else so i could avoid elseif above
        $curVisited{$neighbourI} = 1;
      }
    }

  }
  Log(SLog::ALG_EVAL_GSWN_SUM_VN, "Node: $nodeI Cost: $neighbourCost VisNodes: %visitedNodes");
  # TODO: @visitedNodes not contains %curVisited
  return($neighbourCost, %curVisited);
}


sub preparePreVisited {
  my %preVisited;
  foreach my $node (keys %dG) {
    $preVisited{$node} = 0;
  }
  return %preVisited;
}

sub randIndex {
  # return int(rand((keys %dG) - 1)) + 1;
  # my %di = %_;
  # my %di  = %{shift()};
  my (%di) = @_;
  # our %di; local *di = \$_[0];
  # my %di = %diR;
  # return int(rand(keys %di)) + 1;
  # return $di{(keys %di)[rand keys %di]};
  # my $r = rand keys %di;
  # my @da = keys %di;
  # my $k = $di[rand @da];
  my @hash_keys = keys %di;
  my $k = $hash_keys[rand @hash_keys];
  return $k;
}

sub isResOk {
  my @curResult = @_;
  if (@curResult == 0) {return 0;} # if empty = not executed yet, then result is not ok
  # TODO other ifs ...
  # here i should have some conditions like:
  # - compare max length of path resulting subgraph
  # - compare max tree width of resulting subgraph (if its similar to tree)
  return 1;
}


sub isSufficient {
  my %res = %{shift()};
  $iteration++;
  if ($iteration > $settings{maxIteration}) {
    Log(SLog::ALG_EVO_LOOP_ISSIUF, "Overiterated - ending ($iteration > $settings{maxIteration})");
    return 1;
  }

  if (defined $res{'avgCost'} && defined $bestResSoFar{'avgCost'}) {
    if ($res{'avgCost'} > $bestResSoFar{'avgCost'}) {
      Log(SLog::ALG_EVO_LOOP_ISSIUF, "Continuing ... res.avgCost: $res{'avgCost'} > bestResSoFar.avgCost: $bestResSoFar{'avgCost'}");
      %bestResSoFar = %res;
      return 0;
    }
    else {
      # TODO if 1) X% of current population match previous genome, X>90?
      $sameResCnt++;
      if ($sameResCnt > $settings{maxSameResCnt}) {
        Log(SLog::ALG_EVO_LOOP_ISSIUF, "Breaking... at $iteration with res.avgCost: $res{'avgCost'}");
        Log(SLog::ALG_EVO_LOOP_ISSIUF, \%res);
        return 1;
      }
      return 0;
    }
  }
  else {
    %bestResSoFar = %res;
    return 0;
  }
  # if ($iteration > $dN * $dM) {
  # max iterations ...
  # Log(SLog::ALG_EVO_LOOP_ISSIUF, "Too many iterations");
  # return 1; #
  # }
  Log(SLog::ALG_EVO_LOOP_ISSIUF, "Wow, not sufficient conditioning");
  return 1;
}

sub randomizeGenes {
  my %genes = (
    'startPos'      => randIndex(%dG),                     # alely [1..4]
    'nextPosOffset' => randIndex(%dG) % (scalar(%dG) - 1), # alely [1..3]
    # TODO ^this will not work for more than 2 locations .. '(scalar(%dG) - 1)' should be replaced with freeNodesCnt
    # 'reuseCovered'  => int(rand(1)),                               # TODO
  );
  return %genes;
}

# this subroutine is for test only and I used it hardcoded while debugging...
sub testGenes {
  my %genes = (
    'startPos'      => 4, # alely [1..4]
    'nextPosOffset' => 1, # alely [1..3]
    # 'reuseCovered'  => int(rand(1)),                               # TODO
  );
  return %genes;
}


sub evalEvolutionAlg {
  # my %genes = (
  #   startPos      => 1, # alely [1..4]
  #   nextPosOffset => 1, # alely [1..3]
  # );

  # Population will be always 10
  my $populationMax = 3;
  my $populationCnt = $populationMax;
  my %populationRes = ();

  # generovani(P_0)
  for (1 .. $populationCnt) {
    # my %genes = randomizeGenes();
    my %genes = $settings{isTest} ? testGenes() : randomizeGenes();
    my %singleRes = getCostWithMultiNodes($dV, \%genes);
    $populationRes{$_} = \%singleRes;
  }

  my $t = 0;
  while (!isSufficient(\%populationRes)) {
    # while 'není splněna podmínka ukončení' do
    # vyhodnocení(P_t) # happens in isSufficient();
    # P′_t=výběr(Pt);(strategie výběru)
    # Contest... or smth like that, for simplicity sake I will just take results better then avg ..
    my $totalCost = 0;
    for (1 .. $populationCnt) {
      $totalCost += $populationRes{$_}{'totalCost'};
    }
    my $avgCost = $totalCost / $populationCnt;
    my %desiredParents = ();
    Log(SLog::ALG_EVO_LOOP, "All parents:");
    Log(SLog::ALG_EVO_LOOP, \%populationRes);
    for (1 .. $populationCnt) {
      if ($populationRes{$_}{'totalCost'} >= $avgCost) {
        $desiredParents{$_} = $populationRes{$_};
      }
    }
    Log(SLog::ALG_EVO_LOOP, "AvgCost: $avgCost ; Desired parents:");
    Log(SLog::ALG_EVO_LOOP, \%desiredParents);

    # P′t=reprodukce(P′t);(strategie reprodukce)
    # reproduction, crossbreeding
    my %newPopulationsGenes = ();
    for (1 .. $populationCnt) {
      my $infoMsg = "";
      my $parent1Genes;
      if ($settings{reproductionCrossbreeding}) {
        $infoMsg .= "Crossbreeding genes for $_: \n";
        # TODO dont use same parent..?
        $parent1Genes = $desiredParents{randIndex(%desiredParents)}{'genes'};
        my $parent2Genes = $desiredParents{randIndex(%desiredParents)}{'genes'};

        foreach my $geneName (keys %$parent2Genes) {
          if (int(rand(2)) % 2) {
            $newPopulationsGenes{$_}{$geneName} = $parent1Genes->{$geneName};
          }
          else {
            $newPopulationsGenes{$_}{$geneName} = $parent2Genes->{$geneName};
          }
          $infoMsg .= "  - $geneName: => ";
          $infoMsg .= "$newPopulationsGenes{$_}{$geneName}\n";
        }
      }
      if ($settings{reproductionMutation}) {
        # mutation, why as well as crossbreeding? Because with low initial population I can get really bad results and so
        # never get to the correct result, this brings some possibility of correcting...
        $infoMsg .= "Mutating genes for $_: \n";
        # 'keys %$parent1Genes' will work because number of genes is always the same (not like in real life..)
        foreach my $geneName (keys %$parent1Genes) {
          if (int(rand($settings{mutationRandMax})) >= $settings{mutationRandThreshold}) {
            # pst 1/10
            my $newRandGeneVal = 1;
            if ($geneName eq 'startPos') {
              $newRandGeneVal = randIndex(%dG);
            }
            elsif ($geneName eq 'nextPosOffset') {
              $newRandGeneVal = randIndex(%dG) % (scalar(%dG) - 1);
            }
            $infoMsg .= "  - $geneName: " . $newPopulationsGenes{$_}{$geneName} . " => " . $newRandGeneVal . "\n";
            $newPopulationsGenes{$_}{$geneName} = $newRandGeneVal;
          }
        }
      }
      Log(SLog::ALG_EVO_LOOP, $infoMsg);
    }

    # vyhodnocení(P′t); # there is no comment in lectures :(
    # Is this the same as generating for P_0?
    # I think this is just picking some parents with children, but I am not sure...
    # What I do is replace all parents, but (some) genes persist so results are the same.
    # Pt+1=nahrazení(Pt,P′t); # idk what to do with this one considering the previous one ...
    %populationRes = ();

    for (1 .. $populationCnt) {
      my $genes = $newPopulationsGenes{$_};
      my %singleRes = getCostWithMultiNodes($dV, \%{$genes});
      $populationRes{$_} = \%singleRes;
    }
    $populationRes{'totalCost'} = $totalCost;
    $populationRes{'avgCost'} = $avgCost;

    Log(SLog::ALG_EVO_LOOP, "Iter: " . $iteration . " AvgCost: " . $avgCost . " Res:");
    Log(SLog::ALG_EVO_LOOP, "newPopulationsGenes");
    Log(SLog::ALG_EVO_LOOP, \%newPopulationsGenes);
    Log(SLog::ALG_EVO_LOOP, "populationRes");
    Log(SLog::ALG_EVO_LOOP, \%populationRes);
    # (strategie náhrady)t=t+1
  }
}

# Let's go in to commit the old sin...
main();
#
