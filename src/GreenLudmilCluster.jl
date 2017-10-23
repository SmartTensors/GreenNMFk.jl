# Pairwise Cosine distance
# 	d_st = 1 - (x_s * x_t') / sqrt((x_s * x_s') * (x_t * x_t'))
function pdist_cosine(x, y)
	return (1 - dot(vec(x), vec(y)) / (norm(x) * norm(y)))
end

"""
Ludmil Cluster

$(DocumentFunction.documentfunction(ludmil_cluster;
argtext=Dict("Sources3D"=>"",
			"colSources"=>"",
			"clusterRepeatMax"=>"")))

Returns:
- idx
- idx_r
- Center
"""
function ludmil_cluster(Sources3D, colSources, clusterRepeatMax)

	allProcesses = Sources3D
    
    numberOfPoints = size(allProcesses, 1)
	numberOfProcesses = size(allProcesses, 2)
	globalIter =  size(allProcesses, 3)
	centroids = allProcesses[:, :, 1]

	idx = zeros(numberOfProcesses, globalIter)

	for clusterIt = 1:clusterRepeatMax
		for globalIterID = 1:globalIter

			processesTaken = zeros(numberOfProcesses, 1)
			centroidsTaken = zeros(numberOfProcesses, 1)

			for currentProcessID = 1:numberOfProcesses
				distMatrix = ones(numberOfProcesses, numberOfProcesses) * 100

				for processID = 1:numberOfProcesses
					for centroidID = 1:numberOfProcesses
						if ( (centroidsTaken[centroidID] == 0) && ( processesTaken[processID] == 0) )
							distMatrix[processID, centroidID] = pdist_cosine(allProcesses[:, processID, globalIterID]', centroids[:,centroidID]') #pdist( [ allProcesses[:, processID, globalIterID]' ; centroids[:,centroidID]' ] , 'cosine')
						end
					end
                end

                First = find(distMatrix .== minimum(distMatrix[:]))
				minProcess, minCentroid = ind2sub(size(distMatrix), First[1])
				processesTaken[minProcess] = 1
				centroidsTaken[minCentroid] = 1
				idx[minProcess, globalIterID] = minCentroid
			end
        end

        #idx(:,1) = 1:numberOfProcesses;     %%%%%% HACK
		centroids = zeros(numberOfPoints, numberOfProcesses)
		for centroidID = 1:numberOfProcesses
			for globalIterID = 1:globalIter
				centroids[:, centroidID] = centroids[:, centroidID] + allProcesses[:, find( idx[:, globalIterID] .== centroidID ), globalIterID]
			end
		end
		centroids = centroids ./ globalIter

	end
	
	idx_r = (reshape(idx, numberOfProcesses * globalIter, 1))
    
    Center  = zeros(numberOfProcesses,3)
    
    for d = 1:numberOfProcesses
        Center[d,:] = mean(colSources[find(idx_r == d ),:])
    end

	return idx, idx_r, Center 
end