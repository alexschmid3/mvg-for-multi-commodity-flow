
using Luxor, Colors

#--------------------------------------------------------#

function drawmap(drawingname, mcfinstance, xdim, ydim)

    #Get correct scale
    mult = xdim	 
    shift = xdim/2
        
    #Format and transform latitude and longitude coordinates of each pit stop
    pointDict = Dict()
    listofpoints = []
    listofpoints_labels = []
    for l in 1:numnodes
        xcoord, ycoord = mcfinstance.coordinates[l,1], mcfinstance.coordinates[l,2]
        transformedcoords = (mult*xcoord-shift, mult*ycoord-shift)
        pointDict[l] = Point(transformedcoords)
        push!(listofpoints, transformedcoords)
        push!(listofpoints_labels, [transformedcoords, string(l)])
    end
    locationPoints = Point.(listofpoints)

    #--------------------------------------------------------#

    #Calculate arcs
    arcList = []	
    for a in 1:mcfinstance.numarcs
        i, j = mcfinstance.arcs[a]
        startPoint = locationPoints[i]
        endPoint = locationPoints[j]
        thickness = 3
        push!(arcList, (startPoint, endPoint, (0,0,0), thickness, "solid"))
    end

    #--------------------------------------------------------#

    #Create new drawing
    Drawing(xdim, ydim, drawingname)
    origin()
    background("white")

    #Draw the arcs
    for i in arcList
        #Set arc attributes
        setline(i[4])
        setcolor(i[3])
        setdash(i[5])

        #Draw the arc line
        line(i[1], i[2] , :stroke)
        
        #Calculate the rotation and placement of the arrowhead
        theta = atan((i[2][2] - i[1][2])/(i[2][1] - i[1][1]))
        dist = distance(i[1], i[2])
        arrowhead = (1-8/dist)*i[2] + (8/dist)*i[1] #center of arrowhead positioned 8 pixels from the end node

        #Rotate the arrowhead appropriately
        if i[1][1] >= i[2][1]
            local p = ngon(arrowhead, min(8, i[4]*2), 3, theta - pi , vertices=true)
        else
            local p = ngon(arrowhead, min(8, i[4]*2), 3, theta , vertices=true)
        end

        #Draw the arrowhead
        poly(p, :fill,  close=true)
    end

    #Draw the pit stop nodes
    setcolor("red")
    circle.(locationPoints, 16, :fill)
    setcolor("black")
    setline(3)
    circle.(locationPoints, 16, :stroke)

    #Add pit stop labels
    #for item in listofpoints_labels
    #	label(item[2], :E , Point(item[1]))
    #end

    #--------------------------------------------------------#

    finish()
    preview()

end