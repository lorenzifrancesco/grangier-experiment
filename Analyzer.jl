module Analyzer
using Plots
using Printf
import PyPlot

export analyze, main, loader

pyplot()
default(show = true)
PyPlot.clf()
println(PyPlot.backend)

function loader()
    println("Loading...")
    s = "./tags.txt"
    a = readlines(s)
    for y in a
        filter(x -> !isspace(x), y)
    end
    i=0
    b = Array{Int, 2}(undef, 2, length(a))

    b[1, :] = [parse(Int, split(x, ";")[1]) for x in a]
    b[2, :] = [parse(Int, split(x, ";")[2]) for x in a]


    tags = Array{Int, 2}(undef, 3, length(b))
    fill!(tags, 0)
    println(typeof(tags))
    k = Array{Int, 1}(undef, 3) # k[i] will be the total count of trigger events on channel i
    fill!(k, 1)
    i=0
    cnt = 0
    aft = Array{Int, 1}(undef, 3)
    fill!(aft, 0)
    for i = 1:length(a)
        # if (i>7)
        #     println("condition: ", b[1, i] ," < ", tags[ b[2, i]-1, k[b[2, i]-1] - 1 ], "  ?")
        # end
        if (i<8 || b[1, i] - 3900 > tags[ b[2, i]-1, k[b[2, i]-1] - 1 ])
            tags[ b[2, i]-1, k[b[2, i]-1] ] = b[1, i]
            k[b[2, i] - 1] += 1
        else
            # println("Afterpulse on CH-", b[2, i] - 1)
            aft[b[2, i] - 1] +=1
        end
    end

    println("Number of valid hits")
    @printf("\t n. of transmitted hits   : %6d \n", k[1])
    @printf("\t n. of reflected hits     : %6d \n", k[2])
    @printf("\t n. of gate hits          : %6d \n", k[3])
    println("T+R = ", k[1]+k[2], ", G = ", k[3])
    println("Number of afterpulses:")
    @printf("\t chan 1 - transmitted (2) : %6d \n", aft[1])
    @printf("\t chan 2 - reflected (3)   : %6d \n", aft[2])
    @printf("\t chan 3 - gate (4)        : %6d \n", aft[3])
    println("Percentage of afterpulses")
    @printf("\t chan 1 - transmitted (2) : %4.1f %% \n", aft[1]/k[1] * 100)
    @printf("\t chan 2 - reflected (3)   : %4.1f %% \n", aft[2]/k[2] * 100)
    @printf("\t chan 3 - gate (4)        : %4.1f %% \n", aft[3]/k[3] * 100)
    return (tags, k);
end

function analyze((tags, k))
    println("Analyzing...")
    machine_time = 80.955e-12

    g1 = -1         # BE CAREFUL : NOT A REAL GATE EVENT
    g2 = tags[3, 1]
    diff1 = Array{Int, 1}(undef, k[1])
    diff2 = Array{Int, 1}(undef, k[2])
    fill!(diff1, 0)
    fill!(diff2, 0)
    n = 1

    # Retarded gate method - positive diff
    # for i = 2:k[3]
    #     while (tags[1, n]<g2 && n<k[1])
    #         diff1[n] = g2 - tags[1, n]
    #         n += 1
    #     end
    #     g2 = tags[3, i]
    # end

    # Anticipated gate method - positive diff
    n = 8
    for i = 2:k[3]
        while (tags[1, n]<g2 && n<k[1])
            diff1[n] = tags[1, n] - g1
            n += 1
        end
        g1 = g2
        g2 = tags[3, i]
    end
    diff1 = diff1[8:length(diff1)]
    
    # Minimum distance method
    # for i = 2:k[3]
    #     while (tags[1, n]<g2 && n<k[1])
    #         if ((tags[1, n] - g1) <  (g2 - tags[1, n]))
    #             diff1[n] = tags[1, n] - g1
    #         else
    #             diff1[n] = tags[1, n] - g2
    #         end
    #         n += 1
    #     end
    #     g1 = g2
    #     g2 = tags[3, i]
    # end

    # unreal difference filter

    println("maximum difference CH1     : ", maximum(diff1))
    println("minimum difference CH1     : ", minimum(diff1))
    @printf("maximum time difference CH1 (ns)  : %10.4f \n", maximum(diff1)/machine_time * 1e9)
    @printf("minimum time difference CH1 (ns)  : %10.4f \n", minimum(diff1)/machine_time * 1e9)

    filter!(x-> (x< 25e-9/machine_time), diff1)
    println("Filtered differences are n. : ", length(diff1))
    println("Which is ")

    mod = Int(floor(0.2e-9/machine_time))
    max_diff1 = maximum(diff1)
    min_diff1 = minimum(diff1)
    bin_num = Int(floor((max_diff1-min_diff1) / mod)) + 1
    bias = Int(floor(-min_diff1/mod))
    # sort!(diff1)

    hist = Array{Int, 1}(undef, bin_num)
    println("hist length: ", bin_num)
    
    fill!(hist, 0)
    i = 1
    while (i<=length(diff1))
        hist[Int(floor((diff1[i] - min_diff1) / mod))+1] += 1
        i += 1
    end

    # for i=1:length(hist)
    #     @printf("%20.18f - %d \n", i*80e-12, hist[i])
    # end
    println("maximum frequency at: ", (findmax(hist)[2]-bias)*mod)
    if (length(hist)<600)
        println("Plotting...")
        first_n = length(hist)
        fig = PyPlot.figure()
        fig = Plots.bar((-bias:first_n-bias)*mod*machine_time*1e9,
                         hist[1:first_n],
                         show=true,
                         xlabel = "difference to ** gate event (ns)",
                         ylabel = "Frequency")
        display(fig)
    else
        println("Showing...")
        for i=1:length(hist)
            @printf("%20.10f - %d \n", i, hist[i])
        end
    end

    # COINCIDENCES EVALUATION


    return 0
end

function poisson_pdf(x, mu)
    return (mu)^x/factorial(big(x)) * exp(-mu)
end

function bose_einst_pdf(x, mu)
    return 1/(mu + 1) * (mu/(mu+1))^x
end

function main()
    println("Nothing to do...")
end
end