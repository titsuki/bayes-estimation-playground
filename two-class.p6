use v6;

sub fx($x, $mu = 0, $sigma = 1) {
    1 / sqrt(2 * pi * $sigma) * exp(- ($x - $mu) ** 2 / (2 * $sigma));
}

my Int $num-classes = 2;
my Int $num-samples = 20;
my @data;
my @class-probs[$num-samples;$num-classes];
my @Pc[$num-classes];
my Hash @param[$num-classes];
my @log-likelihood;

@data.push(rand) for ^10;
@data.push(rand + 5.0) for ^10;

for ^$num-classes -> $class {
    @class-probs[$_;$class] = rand for ^@data;
}

for ^@data -> $data-i {
    my $sum = (^$num-classes).map({@class-probs[$data-i;$_]}).sum;
    @class-probs[$data-i;$_] /= $sum for ^$num-classes;
}

for ^1000 -> $iter {
    # M-step
    for ^$num-classes -> $class {
        @Pc[$class] = (^@data).map({@class-probs[$_;$class]}).sum / $num-samples;
    }
    
    for ^$num-classes -> $class {
        my @probs = (^@data).map({@class-probs[$_;$class]});
        @param[$class]<mu> = (@probs Z* @data).sum / @probs.sum;
        @param[$class]<sigma> = (@probs Z* (@data.map({ ($_ - @param[$class]<mu>) ** 2 }))).sum / @probs.elems;
    }
    
    # E-step
    for ^@data -> $data-i {
        for ^$num-classes -> $class {
            @class-probs[$data-i;$class]
            = exp(@Pc[$class].log + fx(@data[$data-i], @param[$class]<mu>, @param[$class]<sigma>).log);
        }
        my $tmp = (^$num-classes).map({ @class-probs[$data-i;$_] }).sum;
        for ^$num-classes -> $class {
            @class-probs[$data-i;$class] /= $tmp;
        }
    }

    # likelihood
    my @logsum[$num-classes];
    for ^$num-classes -> $class {
        @logsum[$class] = (^@data).map({@class-probs[$_;$class]}).sum.log;
    }
    @log-likelihood[$iter] = @logsum.map(*.exp).sum.log;

    @log-likelihood[$iter].say;
    @param.say;
    @Pc.say;
    if $iter > 0 {
        last if @log-likelihood[$iter] + 1e-10 < @log-likelihood[$iter - 1];
    }
}
