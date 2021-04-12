using Distributions

# #-------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------
# # This file contains functions that are used for various climate projection and SC-CO2 calculations.
# #-------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------



#######################################################################################################################
# SAMPLE PAGE2009 UNCERTAIN PARAMETERS
#######################################################################################################################
# Description: This function samples the uncertain parameters specified for PAGE2009 (assuming 
#              independence between different parameters) and stores them in a single dataframe.
#
# Function Arguments:
#
#       n_samples = The number of total parameter samples to take.
#----------------------------------------------------------------------------------------------------------------------

function sample_page_parameters(n_samples::Int)

	# Create a vector of PAGE2009 uncertain parameter names.
	page_param_names = ["stay_fractionCO2emissionsinatm", "res_CO2atmlifetime", "air_CO2fractioninatm", "ccf_CO2feedback", "ccfmax_maxCO2feedback",
				        "tcr_transientresponse", "frt_warminghalflife", "pole_polardifference", "rlo_ratiolandocean",
				        "sltemp_SLtemprise", "sla_SLbaselinerise", "sltau_SLresponsetime", "s0_initialSL",
				        "d_sulphateforcingbase", "ind_slopeSEforcing_indirect"]

	# Initialize dataframe to hold parameter samples.
	param_sample = DataFrame(fill(Float64, length(page_param_names)), Symbol.(page_param_names), n_samples)

	#---------------------------------------------------
	# Get samples of each uncertain PAGE2009 parameter.
	#---------------------------------------------------

	# Fraction of CO₂ emissions remaining in atmosphere.
	param_sample.stay_fractionCO2emissionsinatm = rand(TriangularDist(0.25, 0.35, 0.30), n_samples)
	# Half life of CO₂ atmospheric residence.
	param_sample.res_CO2atmlifetime = rand(TriangularDist(50.0, 100.0, 70.0), n_samples)
	# Percent of CO₂ emitted to air.
	param_sample.air_CO2fractioninatm = rand(TriangularDist(57.0, 67.0, 62.0), n_samples)
	# Stimulation of CO₂ concentration.
	param_sample.ccf_CO2feedback = rand(TriangularDist(4.0, 15.0, 10.0), n_samples)
	# CO₂ stimulation limit.
	param_sample.ccfmax_maxCO2feedback = rand(TriangularDist(30.0, 80.0, 50.0), n_samples)
	# Transient climate response.
	param_sample.tcr_transientresponse = rand(TriangularDist(1.0, 2.8, 1.3), n_samples)
	# Half-life of global warming.
	param_sample.frt_warminghalflife = rand(TriangularDist(10.0, 65.0, 30.0), n_samples)
	# Poles excess temperature change over equator.
	param_sample.pole_polardifference = rand(TriangularDist(1.0, 2.0, 1.5), n_samples)
	# Land excess temperature ratio to ocean.
	param_sample.rlo_ratiolandocean = rand(TriangularDist(1.2, 1.6, 1.4), n_samples)
	# Sea level rise with temperature.
	param_sample.sltemp_SLtemprise = rand(TriangularDist(0.7, 3.0, 1.5), n_samples)
	# Sea level asymptote.
	param_sample.sla_SLbaselinerise = rand(TriangularDist(0.5, 1.5, 1.0), n_samples)
	# Half-life of sea level rise.
	param_sample.sltau_SLresponsetime = rand(TriangularDist(500.0, 1500.0, 1000.0), n_samples)
	# Sea level rise in base year (1849).
	param_sample.s0_initialSL = rand(TriangularDist(-0.04795351, 0.05204649, 0.00204649), n_samples)
    # Sulphate direct (linear) effect in base year (1849).
	#param_sample.d_sulphateforcingbase = rand(TriangularDist(-0.396025801, 0.0, -0.003974199), n_samples)
#	param_sample.d_sulphateforcingbase = rand(TriangularDist(-0.0378, 0.0, -0.003974199), n_samples)
	param_sample.d_sulphateforcingbase = rand(TriangularDist(-0.1, 0.0, -0.003974199), n_samples)
           #!!!! ORIGINAL PAGE VERSION param_sample.d_sulphateforcingbase = rand(TriangularDist(-0.8, -0.2, -0.4), n_samples)
	# Sulphate indirect (log) effect for a doubling.
	       # !!! ORIGINAL PAGE VERSION param_sample.ind_slopeSEforcing_indirect = rand(TriangularDist(-0.8, 0.0, -0.4), n_samples)
	param_sample.ind_slopeSEforcing_indirect = rand(TriangularDist(-0.4, 0.0, -0.060809361), n_samples)

	# Return sampled parameters.
	return param_sample
end



#######################################################################################################################
# SAMPLE FUND UNCERTAIN PARAMETERS
#######################################################################################################################
# Description: This function samples the uncertain parameters specified for FUND (assuming independence between 
#              different parameters) and stores them in a single dataframe.
#
# Function Arguments:
#
#       n_samples = The number of total parameter samples to take.
#----------------------------------------------------------------------------------------------------------------------

function sample_fund_parameters(n_samples::Int)

    # Create a vector of FUND uncertain parameter names.
    fund_param_names = ["terrco2sens", "lifeco2", "lifeco3", "lifeco4", "lifeco5",
                        "lifetempconst", "lifetemplin", "lifetempqd", "climatesensitivity",
                        "lifesea", "seas"]

    # Initialize dataframe to hold parameter samples.
    param_sample = DataFrame(fill(Float64, length(fund_param_names)), Symbol.(fund_param_names), n_samples)

    #---------------------------------------------------
    # Get samples of each uncertain PAGE2009 parameter.
    #---------------------------------------------------

    # Note: Distributions from FUND documentation and Excel tables from previous model versions (when pdfs were't well documented) (https://fundjl.readthedocs.io/en/latest/science.html#atmosphere-and-climate)

    # Parameter in terrestrial biosphere.
    param_sample.terrco2sens = rand(Gamma(4.92255, 662.834125760029), n_samples)
    # Carbon decay time in carbon cycle Box #2.
    param_sample.lifeco2 = rand(Truncated(Normal(363.0, 90.75),0.0,Inf), n_samples)
    # Carbon decay time in carbon cycle Box #3.
    param_sample.lifeco3 = rand(Truncated(Normal(74.0, 18.5),0.0,Inf), n_samples)
    # Carbon decay time in carbon cycle Box #4.
    param_sample.lifeco4 = rand(Truncated(Normal(17.0, 4.25),0.0,Inf), n_samples)
    # Carbon decay time in carbon cycle Box #5.
    param_sample.lifeco5 = rand(Truncated(Normal(2.0, 0.5),0.0,Inf), n_samples)
    # Parameter for temperautre e-folding time.
    param_sample.lifetempconst = rand(Normal(-42.7421952811, 0.12072054186561), n_samples) # Mean is current FUND default, with σ from previous version pdf in documented spreadsheet (only place I could find this info).
    # Parameter for temperautre e-folding time.
    param_sample.lifetemplin= rand(Normal(29.0603120788, 0.0256648976380545), n_samples)# Mean is current FUND default, with σ from previous version pdf in documented spreadsheet (only place I could find this info).
    # Parameter for temperautre e-folding time.
    param_sample.lifetempqd= rand(Normal(0.0014564222822, 0.001568), n_samples)# Mean is current FUND default, with σ from previous version pdf in documented spreadsheet (only place I could find this info).
    # Equilibrium climate sensitivity.
    param_sample.climatesensitivity = rand(Truncated(Gamma(6.47815626, 0.547629469),0.0,Inf), n_samples)
    # Sea level rise e-folding time.
    param_sample.lifesea = rand(TriangularDist(250.0, 1000.0, 500.0), n_samples)
    # Sea level rise sensitivity to temperature.
    param_sample.seas = rand(Gamma(6.0, 0.4), n_samples)

    # Return sampled parameters.
    return param_sample
end



######################################################################################################################
# SAMPLE U.S. CLIMATE SENSITIVTY VALUES
#######################################################################################################################
# Description: This function samples the equilibrium climate sensitivity distribution used for official U.S. social
#              cost of carbon estimates (Roe & Baker distribution with specific parameterization).
#
# Function Arguments:
#
#       n_samples = The number of total parameter samples to take.
#----------------------------------------------------------------------------------------------------------------------

function sample_us_climate_sensitivity(n_samples::Int)

    # Create a distribution representing uncertainty in system feedbacks.
    feedback_dist = Truncated(Normal(0.6198, 0.1841), -0.2, 0.88)

    # Create U.S. climate sensitivity sample.
    ecs_sample = 1.2 ./ (1 .- rand(feedback_dist, n_samples))

    # Return sampled values.
    return ecs_sample
end



#######################################################################################################################
# CREATE CLIMATE PROJECTION CREDIBLE INTERVALS
#######################################################################################################################
# Description: This function calculates the upper and lower credible interval ranges for the climate projections carried
#              out with each posterior parameter sample. It does so for two different % levels.
#
# Function Arguments:
#
#       years          = The years the model projection is carried out for.
#       model_result   = An array of model projections (each row is a new projection, each column is a new year).
#       conf_1_percent = First percentage value to calculate credible interval over (i.e. 0.95 corresponds to 95% credible interval).
#       conf_2_percent = Second percentage value to calculate credible interval over.
#----------------------------------------------------------------------------------------------------------------------

function get_confidence_interval(years, model_result, conf_1_percent, conf_2_percent)

    # Set intervals for quantile function and calculate total number of years in results.
    α1      = 1-conf_1_percent
    α2      = 1-conf_2_percent
    n_years = length(years)

    # Initialize dataframe of results with a column of years, mean results, and missing values.
    ci_results = DataFrame(Year=years, Mean=vec(mean(model_result, dims=1)), Lower1=zeros(Union{Missing,Float64}, n_years), Upper1=zeros(Union{Missing,Float64}, n_years), Lower2=zeros(Union{Missing,Float64}, n_years), Upper2=zeros(Union{Missing,Float64}, n_years))

    # Calculate credible intervals for each year (CI for years with 'missing' values also set to 'missing').
    for i in 1:n_years
        if all(x-> x !== missing, model_result[:,i])
            ci_results.Lower1[i] = quantile(model_result[:,i], α1/2)
            ci_results.Upper1[i] = quantile(model_result[:,i], 1-α1/2)
            ci_results.Lower2[i] = quantile(model_result[:,i], α2/2)
            ci_results.Upper2[i] = quantile(model_result[:,i], 1-α2/2)
        else
            ci_results.Lower1[i] = missing
            ci_results.Upper1[i] = missing
            ci_results.Lower2[i] = missing
            ci_results.Upper2[i] = missing
        end
    end

    # Rename columns to have names specific to user-provided credible interval percentages.
    rename!(ci_results, :Lower1 => Symbol(join(["LowerConf" conf_1_percent], '_')))
    rename!(ci_results, :Upper1 => Symbol(join(["UpperConf" conf_1_percent], '_')))
    rename!(ci_results, :Lower2 => Symbol(join(["LowerConf" conf_2_percent], '_')))
    rename!(ci_results, :Upper2 => Symbol(join(["UpperConf" conf_2_percent], '_')))

    return ci_results
end



#######################################################################################################################
# SIMULATE STATIONARY AR(1) PROCESS WITH TIME VARYING OBSERVATION ERRORS.
#######################################################################################################################
# Description: This function simulates a stationary AR(1) process (given time-varying observation errors supplied with
#              each calibration data set) to superimpose noise onto the climate model projections.
#
# Function Arguments:
#
#       n = Number of time periods (years) the model is being run for.
#       σ = Calibrated standard deviation.
#       ρ = Calibrated autocorrelation coefficient.
#       ϵ = Time-varying observation errors.
#----------------------------------------------------------------------------------------------------------------------

function simulate_ar1_noise(n::Int, σ::Float64, ρ::Float64, ϵ::Array{Float64,1})

    # Define AR(1) stationary process variance.
    σ_process = σ^2/(1-ρ^2)

    # Initialize AR(1) covariance matrix (just for convenience).
    H = abs.(collect(1:n)' .- collect(1:n))

    # Calculate residual covariance matrix (sum of AR(1) process variance and observation error variances).
    # Note: This follows Supplementary Information Equation (10) in Ruckert et al. (2017).
    cov_matrix = σ_process * ρ .^ H + Diagonal(ϵ.^2)

    # Return a mean-zero AR(1) noise sample accounting for time-varying observation error.
    return rand(MvNormal(cov_matrix))
end



#######################################################################################################################
# SIMULATE STATIONARY CAR(1) PROCESS WITH TIME VARYING OBSERVATION ERRORS.
#######################################################################################################################
# Description: This function simulates a stationary CAR(1) process (given time-varying observation errors supplied with
#              each calibration data set) to superimpose noise onto the climate model projections.
#
# Function Arguments:
#
#       n              = Number of time periods (years) the model is being run for.
#       α₀             = Calibrated term describing correlation memory of CAR(1) process.
#       σ²_white_noise = Calibrated continuous white noise process variance term.
#       ϵ              = Time-varying observation errors.
#----------------------------------------------------------------------------------------------------------------------

function simulate_car1_noise(n, α₀, σ²_white_noise, ϵ)

    # Indices for full time horizon.
    indices = collect(1:n)

    # Initialize covariance matrix for irregularly spaced data with relationships decaying exponentially.
    H = exp.(-α₀ .* abs.(indices' .- indices))

    # Define the variance of x(t), a continous stochastic time-series.
    σ² = σ²_white_noise / (2*α₀)

    # Calculate residual covariance matrix (sum of CAR(1) process variance and observation error variances).
    cov_matrix = σ² .* H + Diagonal(ϵ.^2)

    # Return a mean-zero CAR(1) noise sample accounting for time-varying observation error.
    return rand(MvNormal(cov_matrix))
end



#######################################################################################################################
# REPLICATE TIME-VARYING OBSERVATION ERRORS FOR PERIODS WITHOUT COVERAGE.
#######################################################################################################################
# Description: This function creates a time-series of observation errors for the entire model time horizon. For years
#              without observation error estimates, the error remains constant at the average of the ten nearest error
#              values in time.
#
# Function Arguments:
#
#       start_year = The first year to run the climate model.
#       end_year   = The final year to run the climate model.
#       error_data = A vector of time-varying observation errors supplied with each calibration data set.
#----------------------------------------------------------------------------------------------------------------------

function replicate_errors(start_year::Int, end_year::Int, error_data)

    # Initialize model years and measurement error vector.
    model_years = collect(start_year:end_year)
    errors = zeros(length(model_years))

    # Find indices for periods that have observation errors.
    err_indices = findall(x-> !ismissing(x), error_data)

    # Replicate errors for all periods prior to start of observations based on average of 10 nearest errors in time.
    errors[1:(err_indices[1]-1)] .= mean(error_data[err_indices[1:10]])

    # Add errors for periods with observations.
    errors[err_indices[1:end]] = error_data[err_indices[1:end]]

    # Replicate errors for all periods after observation data based on average of 10 nearest errors in time.
    errors[(err_indices[end]+1):end] .= mean(error_data[err_indices[(end-9):end]])

    return errors
end



#######################################################################################################################
# LINEARLY INTERPOLATE MODEL RESULTS TO ANNUAL VALUES
#######################################################################################################################
# Description: This function uses linear interpolation to create an annual time series from the model results.
#
# Function Arguments:
#
#       data    = The non-annual model results to be interpolated
#       spacing = Length of time steps between model output.
#----------------------------------------------------------------------------------------------------------------------

function interpolate_to_annual(data, spacing)

    # Create an interpolation object for the data (assume first and last points are end points, e.g. no interpolation beyond support).
    interp_linear = interpolate(data, BSpline(Linear()))

    # Create points to interpolate for (based on spacing term).
    interp_points = collect(1:(1/spacing):length(data))

    # Carry out interpolation.
    return interp_linear[interp_points]
end


