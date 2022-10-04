def single_gene_sim_transients(N_gen_cc_, N_CG_, state_, CG_positions_,
                       spont_gain_params_, coop_gain_params_, coop_maint_params_):

    
    # N_gen_cc_ = number of cell cycles to simulate for

    # write out state data at every update

    
    
    
    state_time_ = [] # record states at specified timepoints
    state_matrix_ = []
    
    
    dat_time_ = [] # for simulated methylation levels
    dat_n_u_ = []
    dat_n_m_ = []

    # initilise system
    time_ = 0.0
    
    # analyse state
    n_u_, n_h_, n_m_ = count_state(state_,N_CG_)

    # update data lists
    dat_time_.append(time_)
    dat_n_u_.append(n_u_)
    dat_n_m_.append(n_m_)
    
    state_time_.append(0)
    state_matrix_ = np.array(state_)

    
    # calculate initial propensities
    coop_gains_propensities_array_, spont_gains_propensities_vector_, \
            coop_maint_propensities_array_, total_propensities_vector_ = \
                initial_propensities(N_CG_, state_, CG_positions_, spont_gain_params_, 
                                         coop_gain_params_, coop_maint_params_)
    
    
    # churn through the simulation
    # total simulation time = no. of gens * no. of cell cycles per gen.
    while True:
        
        # implement timestep using Gillespies direct method
 
        # find indices that sort propensities into ascending order
        propensities_sorted_indicies_ = np.argsort(total_propensities_vector_)

        cumulative_propensities_ = np.cumsum( np.array(total_propensities_vector_)[propensities_sorted_indicies_] )

        # make array listing site indicies to apply propensities-ordering to so that the site 
        # which changes can be identified
        ordered_site_list_ = np.array( np.arange(N_CG_)[propensities_sorted_indicies_] )

        # get a pair of random numbers
        random_pair_ = np.random.rand(2)

        # calc. time incrememnt delta_t_
        delta_t_ = -np.log(random_pair_[0])*(1/cumulative_propensities_[-1])

        
        # check whether have overshot simulation end and if so break
        if (time_ + delta_t_) > N_gen_cc_:
            # record final state at N_gen_cc_ by duplicating last state before updated
            state_time_.append(N_gen_cc_)
            state_matrix_ = np.vstack((state_matrix_,state_matrix_[-1,:])) 
            dat_time_.append(N_gen_cc_)
            dat_n_u_.append(dat_n_u_[-1])
            dat_n_m_.append(dat_n_m_[-1])
            break
            
        # if haven't overshot simulation end will not have broken loop, so, complete the timestep updates

        # cal. which site will change (specified by change_index_)
        propensities_cut_off_ = random_pair_[1]*cumulative_propensities_[-1]

        change_index_ = ordered_site_list_[ np.searchsorted(cumulative_propensities_, propensities_cut_off_) ]


        # update the time
        time_ += delta_t_

        
        # update state
        if state_[change_index_] == 0:
            state_[change_index_] += 2
        elif state_[change_index_] == 2:
            state_[change_index_] -= 2
        else:
            print('ahdla')
        

        # update the propensiteis
        coop_gains_propensities_array_, spont_gains_propensities_vector_, \
            coop_maint_propensities_array_, total_propensities_vector_ = \
            update_propensities(N_CG_, state_, CG_positions_, spont_gain_params_, 
                                         coop_gain_params_, coop_maint_params_, coop_gains_propensities_array_, 
                            spont_gains_propensities_vector_, coop_maint_propensities_array_, change_index_)
            
        # end of timestep
        
        
        
        # write out state data        
        state_time_.append(time_)
        state_matrix_ = np.vstack((state_matrix_,state_))            
        # end of 'write out state data'
            
        # write out mean properties of state
        # whole gene stats.
        # analyse state
        n_u_, n_h_, n_m_ = count_state(state_,N_CG_)
        # update data lists
        dat_time_.append(time_)
        dat_n_u_.append(n_u_)
        dat_n_m_.append(n_m_)


    
    # package sim_outputs
    # 1. 'transients_state'
    # 2. 'transients_meth_level'

    
    transients_state_ = [dat_time_, dat_n_u_, dat_n_m_]
    transients_meth_level_ = [state_ ,  state_time_, state_matrix_]

    
    sim_outputs_ = [transients_state_, transients_meth_level_]
            
    return sim_outputs_

    #######



