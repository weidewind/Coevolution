 use Bio::Phylo::IO;
 use Bio::Phylo::Forest::NodeRole;
 use Bio::Phylo::Forest::TreeRole;
 use Bio::Phylo::Taxa::Taxon;
 use Bio::Phylo::Factory;
 use Bio::Phylo::Listable;

sub test{
print "my test node";
   my $testnode  = $tree -> get_by_name("5_0.600000");
   print "test grand parent ".$testnode ->get_parent ->get_parent->get_name."\n";
   my $node;
   while (1){
    $testnode->visit_depth_first( -pre => sub{ 
    	
    	  $node = shift;
    	 if ($node eq $testnode ->get_parent->get_parent || $node ->is_root){
    	 	print "no!";
    	 	last;
    	 }
   	 print $node ->get_name."\n";
    },
    -pre_sister     => sub {
    	 my $node2 = $_[1]; #current_sister
    	 print "node2 ".$node2->get_name;
    	 my $node1 = ${$testnode->get_sisters}[1]; #that damn sister
    	 print " node1 ".$node1->get_name."\n";
    	 if ($node2->get_name eq $node1->get_name){
    	 	last;
    	 }
    },
    -with_relatives => 1);
   }
 
 }
 