//ccl: Function to print the genotype to logfile
void printGeno(
    wordList genotype,
    label  &i
){
    if (genotype[i] == "I1") {
        i++; Info << "I1";
    } else if (genotype[i] == "I2") {
        i++; Info << "I2";
    } else if (genotype[i] == "*") {
        i++;
        Info << "(";
        printGeno(genotype, i);
        Info << " * ";
        printGeno(genotype, i);
        Info << ")";
    } else if (genotype[i] == "-") {
        i++;
        Info << "(";
        printGeno(genotype, i);
        Info << " - ";
        printGeno(genotype, i);
        Info << ")";
    } else if (genotype[i] == "+") {
        i++;
        Info << "(";
        printGeno(genotype, i);
        Info << " + ";
        printGeno(genotype, i);
        Info << ")";
    } else {
        i++;
        Info << genotype[i-1];
    }

}

//ccl: Function to evaluate the genotype 
volScalarField evalGeno (
    wordList genotype,
    label  &i,
    volScalarField I1,
    volScalarField I2
){
    volScalarField one("one", (I1+1e-10)/(I1+1e-10)); // ugly hack, cast to volScalarField :(
    if (genotype[i] == "I1") {
        i++; 
        return I1;
    } else if (genotype[i] == "I2") {
        i++; return I2;
    } else if (genotype[i] == "*") {
        i++;
        return evalGeno(genotype, i, I1, I2) * evalGeno(genotype, i, I1, I2);
    } else if (genotype[i] == "-") {
        i++;
        volScalarField tmp0("tmp0", evalGeno(genotype, i, I1, I2));
        volScalarField tmp1("tmp1", evalGeno(genotype, i, I1, I2));
        volScalarField tmp("tmp", tmp0 - tmp1);
        return tmp;
    } else if (genotype[i] == "+") {
        i++;
        return evalGeno(genotype, i, I1, I2) + evalGeno(genotype, i, I1, I2);
    } else {
        i++;
        char* pEnd;
        scalar numValue = strtod(genotype[i-1].c_str(),&pEnd);
        return numValue * one;
    }
    return one; // dummy return
}
