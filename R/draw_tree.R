
#this is the function which draws the tree structure

draw_tree = function(out) {

#requires library(ape)

library(ape)
convert_to_newick = function(n) {
	if (n$leaf){
		return(n$label);
	}else{
		return(paste("(",convert_to_newick(n$left),",",convert_to_newick(n$right),")",sep=""));
	}
}

b_tree = function(out) {
	Z = out$Z
	nodes = list();
	Nnodes = dim(Z)[1];
	for(i in 1:leafs) {
		n = list();
		n$label = i;
		n$length = 0.05;
		n$leaf = TRUE;
		nodes[[i]]=	n;
	}
	#build tree
	for(i in 1:Nnodes) {
		#merge according to Z
		left = Z[i,1];
		right = Z[i,2];
		#new node
		n = list();
		n$left = nodes[[left]];
		n$right = nodes[[right]];
		n$leaf = FALSE;
		n$label = i+leafs;
		nodes[[i+leafs]]=n;
	}
	#return nodes
	return(nodes);
}

preorder_labels = function(node) {
	if (node$leaf) {
		return(vector());
	}
	else {
		return(c(node$label,preorder_labels(node$left),preorder_labels(node$right)));
	}
}	
	leafs = dim(out$all_nodes)[2];
	nodes = b_tree(out);
	print(nodes)
	root = nodes[[length(nodes)]];
	stri = paste(convert_to_newick(root),";",sep="");
	print(stri)
	t = read.tree(text=stri);
	plot(t);
	preorder = preorder_labels(root);
	for(i in 1:length(preorder)) {
		nodelabels(preorder[i],i+leafs);
	}
}
