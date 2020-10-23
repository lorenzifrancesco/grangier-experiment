module Analyzer
using Plots
using Printf
import Plotly
import Statistics
export delay_estimator, main, loader, difference_info, gated_counter

Plots.plotly()
default(show = true)
# PyPlot.clf()
# println(PyPlot.backend)
const machine_time = 80.955e-12

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

function delay_estimator((tags, k); mode = "gate_first")
    println("Analyzing...")
    machine_time = 80.955e-12
    diff1 = Array{Int, 1}(undef, k[1])
    diff2 = Array{Int, 1}(undef, k[2])
    fill!(diff1, 0)
    fill!(diff2, 0)
    if mode == "gate_last"
        g1 = -1         # BE CAREFUL : NOT A REAL GATE EVENT
        g2 = tags[3, 1]
        n = 1
        # Retarded gate method - positive diff
        for i = 2:k[3]
            while (tags[1, n]<g2 && n<k[1])
                diff1[n] = g2 - tags[1, n]
                n += 1
            end
            g2 = tags[3, i]
        end

        g1 = -1         # BE CAREFUL : NOT A REAL GATE EVENT
        g2 = tags[3, 1]
        n = 1
        for i = 2:k[3]
            while (tags[2, n]<g2 && n<k[2])
                diff2[n] = g2 - tags[2, n]
                n += 1
            end
            g2 = tags[3, i]
        end
    elseif mode == "gate_first"
        # Anticipated gate method - positive diff
        g1 = -1         # BE CAREFUL : NOT A REAL GATE EVENT
        g2 = tags[3, 1]
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

        g1 = -1         # BE CAREFUL : NOT A REAL GATE EVENT
        g2 = tags[3, 1]
        n = 8
        for i = 2:k[3]
            while (tags[2, n]<g2 && n<k[1])
                diff2[n] = tags[2, n] - g1
                n += 1
            end
            g1 = g2
            g2 = tags[3, i]
        end
        diff2 = diff2[8:length(diff2)]
    else
        # Minimum distance method
        g1 = -100000000         # BE CAREFUL : NOT A REAL GATE EVENT
        g2 = tags[3, 1]
        n = 1
        for i = 2:k[3]
            while (tags[1, n]<g2 && n<k[1])
                if ((tags[1, n] - g1) <  (g2 - tags[1, n]))
                    diff1[n] = tags[1, n] - g1
                else
                    diff1[n] = tags[1, n] - g2
                end
                n += 1
            end
            g1 = g2
            g2 = tags[3, i]
        end
        g1 = -100000000         # BE CAREFUL : NOT A REAL GATE EVENT
        g2 = tags[3, 1]
        n = 1
        for i = 2:k[3]
            while (tags[2, n]<g2 && n<k[1])
                if ((tags[2, n] - g1) <  (g2 - tags[2, n]))
                    diff2[n] = tags[2, n] - g1
                else
                    diff2[n] = tags[2, n] - g2
                end
                n += 1
            end
            g1 = g2
            g2 = tags[3, i]
        end
    end

    max_delay = 10 # [ns]
    max_clicks = max_delay * 1e-9/machine_time
    @printf("PRE-filtering at max delay = %d ns \n ", max_delay)
    # unreal difference filter
    filter!(x-> (x< max_clicks), diff1)
    filter!(x-> (x< max_clicks), diff2)

    difference_info(diff1, diff2, k)
    μ1 = Statistics.mean(diff1)
    μ2 = Statistics.mean(diff2)
    σ1 = sqrt(Statistics.var(diff1 .- μ1))
    σ2 = sqrt(Statistics.var(diff2 .- μ2)) 

    return [μ1, μ2, σ1, σ2]
end

function difference_info(diff1, diff2, k)
    machine_time = 80.955e-12
    println("Difference Info...")
    max_diff1 = maximum(diff1)
    min_diff1 = minimum(diff1)
    max_diff2 = maximum(diff2)
    min_diff2 = minimum(diff2)
    @printf("1) maximum difference      : %10d \n", max_diff1)
    @printf("1) minimum difference      : %10d \n", min_diff1)
    @printf("1) maximum time difference  (ns)  : %10.4f \n", max_diff1 * machine_time * 1e9)
    @printf("1) minimum time difference  (ns)  : %10.4f \n", min_diff1 *machine_time * 1e9)

    @printf("2) maximum difference      : %10d \n", max_diff2)
    @printf("2) minimum difference      : %10d \n", min_diff2)
    @printf("2) maximum time difference  (ns)  : %10.4f \n", max_diff2*machine_time * 1e9)
    @printf("2) minimum time difference  (ns)  : %10.4f \n\n", min_diff2*machine_time * 1e9)

    @printf("1) Percentage of accepted hits : %d / %d = %4.2f \n", length(diff1), k[1], length(diff1)/k[1])
    @printf("2) Percentage of accepted hits : %d / %d = %4.2f\n", length(diff2), k[2], length(diff2)/k[2])

    # Want to show exactly 100 bins in histogram
    mod = Int(ceil(maximum([length(diff1), length(diff2)]) / 10000000))
    mod = 1
    println("Total differences 1 : ", length(diff1))
    println("Total differences 2 : ", length(diff2))

    x_delays1 = (min_diff1:max_diff1)*mod*machine_time*1e9
    x_delays2 = (min_diff2:max_diff2)*mod*machine_time*1e9

    bin_num1 = Int(floor((max_diff1-min_diff1) / mod)) + 1
    println("bins 1: ", bin_num1)
    bias1 = Int(floor(-min_diff1/mod))
    hist1 = Array{Int, 1}(undef, bin_num1)
    fill!(hist1, 0)
    i = 1
    while (i<=length(diff1))
        hist1[Int(floor((diff1[i] - min_diff1) / mod))+1] += 1
        i += 1
    end

    bin_num2 = Int(floor((max_diff2-min_diff2) / mod)) + 1
    bias2 = Int(floor(-min_diff2/mod))
    println("bins 2: ", bin_num2)
    hist2 = Array{Int, 1}(undef, bin_num2)
    fill!(hist2, 0)
    i = 1
    while (i<=length(diff2))
        hist2[Int(floor((diff2[i] - min_diff2) / mod))+1] += 1
        i += 1
    end

    if (length(hist1)<600 && length(hist2)<600)
        println("Plotting...")
        # fig = Plotly.figure()
        fig = Plots.bar(x_delays1,
                         hist1,
                         show=true,
                         xlabel = "difference to ** gate event (ns)",
                         ylabel = "Frequency")
        Plots.bar!(x_delays2, hist2)
        display(fig)
    else
        println("Too long to plot...")
        # for i=1:length(hist1)
        #     @printf("%20.10f - %d \n", i, hist1[i])
        # end
        # for i=1:length(hist2)
        #     @printf("%20.10f - %d \n", i, hist2[i])
        # end
    end
end

# need to decide what method to use
function gated_counter(tags, μ1, σ1, μ2, σ2)
    println("Gated counting...")
    tr_hits = 0
    ref_hits = 0
    coincidences = 0
    n_σ = 3
    for i=1:k[3]
        if (-n_σ*σ1 < (tags[3, i] - tags[1, j] - μ1) < n_σ*σ1)
        end
        if (-n_σ*σ2 < (tags[3, i] - tags[2, k] - μ2) < n_σ*σ2)
        end
    end
end

function main()
    println("Nothing to do...")
end
end