module Analyzer
using Plots
using Printf
import Plotly
import PGFPlots
import Statistics
import ProgressMeter

export delay_estimator, loader, difference_info, gated_counter, single_chan_stat, config

default(show = true)
# PyPlot.clf()
# println(PyPlot.backend)
const machine_time = 80.955e-12

function loader(;aft_filter = true)
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
    if (aft_filter)
        aft_const = 3900
    else
        aft_const = 0
    end
    for i = 1:length(a)
         if (i<8 || tags[ b[2, i]-1, k[b[2, i]-1] - 1 ] + aft_const < b[1, i] )
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

    # max_delay = 7.5 # [ns]
    # max_clicks = max_delay * 1e-9/machine_time
    max_clicks = 100
    max_delay = max_clicks * machine_time / 1e-9
    @printf("PRE-filtering at max delay = %d ns \n ", max_delay)
    # unreal difference filter
    filter!(x-> (x< max_clicks), diff1)
    filter!(x-> (x< max_clicks), diff2)

    filter!(x -> (x>0), diff2)
    # filter!(x -> (x>0), diff1)

    difference_info(diff1, diff2, k)
    μ1 = Statistics.mean(diff1)
    μ2 = Statistics.mean(diff2)
    σ1 = sqrt(Statistics.var(diff1 .- μ1))
    σ2 = sqrt(Statistics.var(diff2 .- μ2)) 

    return [μ1, σ1, μ2, σ2]
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

    @printf("1) Fraction of accepted hits : %d / %d = %4.2f \n", length(diff1), k[1], length(diff1)/k[1])
    @printf("2) Fraction of accepted hits : %d / %d = %4.2f\n", length(diff2), k[2], length(diff2)/k[2])

    # Want to show exactly 100 bins in histogram
    mod = Int(ceil(maximum([length(diff1), length(diff2)]) / 1e4)) # TO BE MODIFIED

    # plot clicks 
    x_delays1 = (min_diff1:mod:max_diff1)
    x_delays2 = (min_diff2:mod:max_diff2)

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
    μ1 = Statistics.mean(diff1)
    μ2 = Statistics.mean(diff2)
    σ1 = sqrt(Statistics.var(diff1 .- μ1))
    σ2 = sqrt(Statistics.var(diff2 .- μ2)) 

    if (length(hist1)<600 && length(hist2)<600)
        println("Plotting...")
        # fig = Plotly.figure()
        n_σ = 2
        fig = Plots.bar(x_delays1,
                         hist1,
                         show=true,
                         xlabel = "absolute difference from gate event (MTU)",
                         ylabel = "frequency", 
                         label = "T", 
                         size = (600, 400))
        Plots.bar!(x_delays2, hist2, label = "R")
        rectangle(w, h, x, y) = Plots.Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

        recr = rectangle(2*n_σ*σ1, maximum([maximum(hist1), maximum(hist2)]), μ1-n_σ*σ1, 0)
        rect = rectangle(2*n_σ*σ2, maximum([maximum(hist1), maximum(hist2)]), μ2-n_σ*σ2, 0)
        # Plots.plot!(recr, linewidth = 2, opacity = 0.1, color=:blue, label=nothing)
        # Plots.plot!(rect, linewidth = 2, opacity = 0.1, color=:red, label=nothing)

        display(fig)
        savefig("./images/delays.pdf")
    else
        println("Too long to plot...")
    end
end

# need to decide what method to use -> we use  GATE -> REFLECTED -> TRANSMITTED
function gated_counter((tags, k), params; mode = "confidence")
    println("Gated counting...")
    μ1 = params[1]
    σ1 = params[2]
    μ2 = params[3]
    σ2 = params[4]

    @printf("mean tramsmitted   : %6.4f \n", params[1])
    @printf("stdd tramsmitted   : %6.4f \n", params[2])
    @printf("mean     reflected : %6.4f \n", params[3])
    @printf("stdd     reflected : %6.4f \n", params[4])
    # N_1 = length(tags[3, :])
    N_1 = 0
    intervals = [2]
    # Gate function (not counting with multiple hits)
    for n_σ in intervals
        max_clicks = 100
        x = 1
        r_hit = false
        refl = 0
        multiple_refl = 0
        y = 1
        t_hit = false
        tran = 0
        multiple_tran = 0
        coincidences = 0
        ## BE CAREFUL INVERTED DATA SET
        if (mode == "confidence") 
            for i=1:length(tags[3, :])-1
                r_hit = false
                t_hit = false  
                while  tags[1, x] < -n_σ*σ1 + tags[3, i] + μ1
                    x += 1
                end
                while -n_σ*σ1 + tags[3, i] + μ1 <= tags[1, x] < +n_σ*σ1 + tags[3, i] + μ1 && tags[1, x] < tags[3, i+1] 
                    r_hit = true
                    x += 1
                end
                if r_hit
                    refl += 1
                end

                while  tags[2, y] < -n_σ*σ2 + tags[3, i] + μ2 
                    y += 1
                end
                while -n_σ*σ2 + tags[3, i] + μ2 <= tags[2, y] < +n_σ*σ2 + tags[3, i] + μ2  && tags[2, y] < tags[3, i+1]
                    t_hit = true
                    y += 1
                end
                if t_hit
                    tran += 1
                end
                if r_hit && t_hit
                    coincidences += 1
                end
                if r_hit || t_hit
                    N_1 += 1
                end
            end
        else
            for i=1:length(tags[3, :])-1
                r_hit = false
                t_hit = false 
                while tags[1, x] < tags[3, i]
                    x += 1
                end
                while tags[3, i] <= tags[1, x] < tags[3, i] + max_clicks
                    r_hit = true
                    x += 1
                end
                if r_hit
                    refl += 1
                end

                while tags[2, y] < tags[3, i]
                    y += 1
                end
                while tags[3, i] <= tags[2, y] < tags[3, i] + max_clicks
                    t_hit = true
                    y += 1
                end
                if t_hit
                    tran += 1
                end
                if r_hit && t_hit
                    coincidences += 1
                end
                if r_hit || t_hit
                    N_1 += 1
                end
            end
        end
        @printf("Measurement with ± σ confidence \n")
        println("sigma = ", n_σ)
        prob_refl = refl / N_1
        prob_tran = tran / N_1
        prob_triple = coincidences / N_1
        α = prob_triple/ (prob_refl * prob_tran)
        @printf("\t gate         hits :  %9d \n", N_1)
        @printf("\t reflected    hits :  %9d \n", refl)
        @printf("\t transmitted  hits :  %9d \n", tran)
        @printf("\t coincidences hits :  %9d \n", coincidences)
        @printf(" ----------------------\n")
        @printf("\t P[double]         : %9.8f \n", prob_refl + prob_tran - 2 *prob_triple)
        @printf("\t P[triple]         : %9.8f \n", prob_triple)
        @printf("\t Alpha             : %9.8f \n", α)

        sigma_r = sqrt(prob_refl*(1-prob_refl)/N_1)
        sigma_t = sqrt(prob_tran*(1-prob_tran)/N_1)
        sigma_c = sqrt(prob_triple*(1-prob_triple)/N_1)
        @printf("p_r variance: %9.8f \n", sigma_r)
        @printf("p_t variance: %9.8f \n", sigma_t)
        @printf("p_c variance: %9.8f \n", sigma_c)
        @printf(" variance: %9.8f \n", sigma_c/(prob_refl*prob_tran) + 
                sigma_r * prob_triple/(prob_refl^2*prob_tran) +
                sigma_t * prob_triple/(prob_refl*prob_tran^2))
    end
end

function config()
    Plots.plotly()
    Plots.default(size=(600, 400), 
    guidefont=("times", 14), 
    legendfont=("times", 14),
    tickfont=("times", 14)
    )
end

function single_chan_stat((tags, k))
    machine_time = 80.955e-12
    bin_num = 1000
    hist = Array{Int, 2}(undef, 3, bin_num)
    fill!(hist, 0)
    bin_step = Array{Int}(undef, 3)
    diff = Array{Int, 2}(undef, 3, maximum(k)-1)
    fill!(diff, 0)
    println(length(tags[3, :]), k)
    maxx = 0
    for chan in [1, 2, 3]
        series = tags[chan, :]
        for i = 1:k[chan]-1
            diff[chan, i] = series[i+1] - series[i]
        end
        filter!(z -> (z>0), diff[chan, :])
        max_diff = maximum(diff[chan, :])
        if max_diff>maxx
            maxx = max_diff
        end
    end
    x_axis = 0:bin_num:maxx
    bin_size = maxx/bin_num
    i = 1
    for chan = [1, 2, 3]
        filter!(z -> (z>0), diff[chan, :])
        for i = 1:k[chan]-2
            hist[chan, Int(ceil(diff[chan, i]/bin_size))] += 1
        end
    end

    fig = Plots.plot((0:bin_num-1)*bin_size,
                     [log10(x) for x in hist[1, :]] ,
                     label = string("trasmitted channel"),
                     show=true,
                     xlabel = "Interarrival time (MTU)",
                     ylabel = "Frequency (log)",
                     size = (600, 400))

    Plots.plot!((0:bin_num-1)*bin_size,
                [log10(x) for x in hist[2, :]],
                label = string("reflected channel"))

    Plots.plot!((0:bin_num-1)*bin_size,
                [log10(x) for x in hist[3, :]],
                label = string("gate channel"))
    @printf("Sum 1 : %5.4f \n", sum(hist[1, :]/sum(hist[1, :])))
    @printf("Sum 2 : %5.4f \n", sum(hist[2, :]/sum(hist[2, :])))
    @printf("Sum 3 : %5.4f \n", sum(hist[3, :]/sum(hist[3, :])))

    savefig(string("./images/single_chan.pdf"))
end

function poisson_moments(mu)
    return [mu, mu, 1/sqrt(mu), 1/mu]
end

function bose_ein_moments(mu)
    sigma = sqrt(mu + mu^2)
    return [mu,
            sigma^2,
            (mu + 3*mu^2 + 2*mu^3)/sigma^3,
            (mu + 10*mu^2 + 18*mu^3 + 9*mu^4)/sigma^4]
end
end