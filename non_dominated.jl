
function add_to_lst(lst, elmt)
    if length(lst) > 0
        insert_idx = -1
        to_remove_idx = []
        if elmt[1] <= lst[1][1] && elmt[2] < lst[1][2]
            # avant le premier
            insert_idx = 1
            push!(to_remove_idx, 1)
        end
        if elmt[1] > lst[length(lst)][1] && elmt[2] < lst[length(lst)][2]
            # après le dernier
            push!(lst, elmt)
        else
            for i = 2:length(lst)
                if insert_idx != -1 && elmt[2] <= lst[i][2]
                    # elmt va être inséré avant, donc on supprime tout élément avec z2 supérieur
                    push!(to_remove_idx, i)
                elseif insert_idx != -1
                    # les éléments suivant on z1 plus grand mais z2 plus petit on peut s'arrêter
                    break
                elseif lst[i][1] < elmt[1] && lst[i][2] < elmt[2]
                    # dominé par un point de la liste, on s'arrête
                    return
                elseif lst[i][1] > elmt[1]
                    insert_idx = i
                    if elmt[2] <= lst[i][2]
                        push!(to_remove_idx, i)
                    end
                elseif lst[i][1] == elmt[1]
                    if elmt[2] < lst[i][2]
                        insert_idx = i
                        push!(to_remove_idx, i)
                    else
                        return
                    end
                end
            end
            if insert_idx != -1
                insert!(lst, insert_idx, elmt)
            end
            for idx_remove in reverse(to_remove_idx)
                if idx_remove >= insert_idx && insert_idx != -1
                    popat!(lst, idx_remove + 1)
                else
                    popat!(lst, idx_remove)
                end
            end
        end
    else
        push!(lst, elmt)
    end
end