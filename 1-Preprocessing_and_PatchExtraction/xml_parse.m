function[xy] = xml_parse(xml_file)
    xDoc = xmlread(xml_file); %get xml doc
    Regions=xDoc.getElementsByTagName('Region'); % get a list of all the region tags

    %get region labels
    for regioni = 0:Regions.getLength-1 %region tags start at 0
        Region=Regions.item(regioni);  %for each region tag
        label = Region.getAttribute('Text'); %ground truth label
        verticies=Region.getElementsByTagName('Vertex'); %get a list of all the vertexes (which are in order)
        xy{1,regioni+1}=zeros(verticies.getLength-1,2); %allocate space for them
        xy{2,regioni+1}=char(label); %grab label for verticies 
        for vertexi = 0:verticies.getLength-1 %iterate through all verticies
            x=str2double(verticies.item(vertexi).getAttribute('X')); %get the x value of that vertex
            y=str2double(verticies.item(vertexi).getAttribute('Y')); %get the y value of that vertex
            xy{1,regioni+1}(vertexi+1,:)=[x,y]; % finally save them into the array
        end   
    end
end


