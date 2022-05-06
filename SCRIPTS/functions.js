function finder(aff,unaff) {
    for(i=0;i<unaff.length;i++){
        if((unaff[i].alts != 0)) { return false; }
    }
    for(i=0;i<aff.length;i++){
        if((aff[i].alts != 1 && aff[i].alts != 2)) { return false; }
    }
	return true
}
