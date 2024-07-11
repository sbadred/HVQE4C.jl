function CheckCommutivity(a::String, b::String, QubitNumber::Int)
result = 1
for i in 1:QubitNumber
    if (a[i] == 'X' && b[i] == 'Y') || (a[i] == 'X' && b[i] == 'Z') ||
       (a[i] == 'Y' && b[i] == 'X') || (a[i] == 'Y' && b[i] == 'Z') ||
       (a[i] == 'Z' && b[i] == 'X') || (a[i] == 'Z' && b[i] == 'Y')
        result *= -1
    end
end
return result == 1
end


 
function CheckInseprability(Pool)
    Set=[Pool[1]]
    i=0
    while i< length(Set)
        for PoolElement in Pool
            if PoolElement ∉ Set
                if !CheckCommutivity(PoolElement,Set[i+1],QubitNumber)
                    push!(Set,PoolElement)
                end
            end 
        end
        i+=1
    end
    if length(Set)< length(Pool)
        return false
    else
        return true
    end
end


function MultiplyStrings(a::String, b::String, QubitNumber::Int)
    result = ""
    for i in 1:QubitNumber
        if a[i] == 'X' && b[i] == 'Y'
            result *= 'Z'
        elseif a[i] == 'X' && b[i] == 'Z'
            result *= 'Y'
        elseif a[i] == 'Y' && b[i] == 'X'
            result *= 'Z'
        elseif a[i] == 'Y' && b[i] == 'Z'
            result *= 'X'
        elseif a[i] == 'Z' && b[i] == 'X'
            result *= 'Y'
        elseif a[i] == 'Z' && b[i] == 'Y'
            result *= 'X'
        elseif a[i] == 'Z' && b[i] == 'Z'
            result *= 'I'
        elseif a[i] == 'X' && b[i] == 'X'
            result *= 'I'
        elseif a[i] == 'Y' && b[i] == 'Y'
            result *= 'I'
        elseif a[i] == 'I'
            result *= b[i]
        elseif b[i] == 'I' && a[i] != 'I'
            result *= a[i]
        end
    end
    return result
end


function GenerateGroup!(starters, QubitNumber)
    tmpSet=[]
    push!(tmpSet,'I'^QubitNumber)
    helpSet=[]
    count=0
    for Pauli in starters
        if Pauli ∈ tmpSet
            count+=1
        end
        if Pauli ∉ tmpSet
            for el in tmpSet
                tmp=MultiplyStrings(Pauli,el,QubitNumber)
                push!(helpSet,tmp)
            end
            #@show tmpSet, helpSet
            L=length(tmpSet)+length(helpSet)
            tmpSet =union(tmpSet, helpSet)
            # if length(tmpSet) < L
            #     break
            # end
            #push!(tmpSet,helpSet)
            helpSet=[]
        end
    end
    @show length(tmpSet)
    return tmpSet, length(tmpSet)
end
 
function GeneratePauliStrings!(set,QubitNumber)
    tmpSet=copy(set)
    MySet=copy(set)
    helpSet=[]

    while length(tmpSet) ≠ 0
        for e1 in tmpSet
            for e2 in MySet
                if !CheckCommutivity(e1,e2,QubitNumber)
                    tmp=MultiplyStrings(e1,e2,QubitNumber)
                    if tmp ∉ MySet
                        ûsh!(helpSet,tmp)
                    end
                end
            end
        end
        MySet=union(MySet, helpSet)
    end
    return MySet
end



function main()
    QubitNumber=8
    starters = ["XZZYYIZY", "YZXIYZYI", "YZIYXZZY", "YIIXYZZY", "YYXXYYXY", "ZXXYIYXY", "ZIIZYYXY", "YYXZZZZY", "ZZZXIIZY", "ZZIXXXYI", "YXIYIZYI", "YXZYZZYI","XZYIZYZY","IIXIIIYI", "YIYIIXZY", "XYIYYYZY", "ZZXIZZYI", "IIIIYYXY","XIXYYIXY"]
    starterList=["XZZYYIZY"]
    
    starters = ["YIYIXZYI", "ZXXIIXYI", "ZIXYXXII", "ZYYIZXYI", "IXXYXIII", "XXXYZZII", "IXXXIXXY", "XXYXXYXY", "YXYXYYXY", "ZXYXXZII", "XYXYXXXY", "YIXZZYZY", "YIZIIYXY", "XIXYYIXY", "XXYZYXYI"]
    starterList =["YIYIXZYI" ]

    starters = [
        "YZXZYZYI", "IZYYIIXY", "YZIXIYYI", "XZYZYIYI", "ZZXXZZXY", "XZYIYIYI", 
        "ZYXZYZZY", "YIIYXIZY", "YIZYIXYI", "YZYZXZYI", "XYZIYYII", "YIXIYIYI", 
        "XIXZXIYI", "ZYYIXIZY", "XXZZYXII", "XZZYIYYI", "ZIYYZIXY", "XIIYYZZY", 
        "XZXIXZYI", "XZZXZXYI", "XIZYYIZY", "IYYIXZZY", "YIIXZYYI", "ZXXIIXYI", 
        "IXIXZXZY", "ZXYIIYYI", "YZIXZYYI", "YZXIYZYI", "XXZIIZXY", "XXIZXYII", 
        "XIZZZXXY", "ZZZXXXYI", "IYIXYIYI", "XIYYXZXY", "XXYXYXXY", "XYXZXYYI", 
        "YYYXZIII"
    ]
    
    
    starterList =["YZXZYZYI"]
    
    extra=[]
    ~,L= GenerateGroup!(starterList, QubitNumber)
    Len=[L]

    for i=2:length(starters)
        push!(starterList,starters[i])
        ~,L= GenerateGroup!(starterList, QubitNumber)
        @show i
        if Len[end] == L
            push!(extra,starters[i])
            deleteat!(starterList, findall(x->x==starters[i],starterList))
        else
            push!(Len,L)
        end
    end

    Group,L= GenerateGroup!(starterList, QubitNumber)
    for i=1:length(Group)
        @show CheckInseprability([Group[i];extra[1]])
    end

    QubitNumber=10
    starters= [
        "YZZYZXYZII", "ZIXYIYIZYI", "ZIYYIZIIXY", "ZYIYIZIXZY", "IZIYYZZXYI", "IIZXZYYIYI", 
        "XIYXZYXXXY", "YZIZYZXXXY", "IYXXXZYXXY", "YYXZXXIXXY", "XZXZXYZYYI", "YZZZZYYXII", 
        "XIYZIYIYII", "IXZZYIYYII", "ZYZIYYYXZY", "XYYYYIZZYI", "YYXIZIIYII", "YIXYIZIIZY", 
        "YIXXXIXXII", "XXYXIYYYYI", "YIYZYYIXYI", "YYXYXIXYZY", "XZYYZIYYYI", "YXXZIIYIXY", 
        "YXYXZXIZZY", "ZYYYXZZZII", "IZXIIXXZZY", "XIIIZZYYZY", "YIYXXIIZXY", "YZZIZYXYII", 
        "XZYYZIZIZY", "YIXZXYIXYI", "ZZIZYYZZXY", "YZYIIXXIXY", "XXXYZYIIZY", "YZIIIYZZXY", 
        "IXYZZYIXXY", "IZIIZIXXXY", "YXYYYIZZYI", "YIXXXIYYII", "IIIXXZIXYI", "XXXYIYZZZY", 
        "YZXYYYIZYI", "IZXIYZIYZY", "XIZYXYIXZY", "XZYIIYIYII", "ZZXXIXXXZY", "YXYYZZIZII", 
        "XZYYIZIZZY", "XIZIYYZIZY", "ZXIYZIYZYI", "XXXXZXXXYI", "ZIXIZIYIII", "YYXIXXZXXY", 
        "IIIXZYYIYI", "YIZIZYZIXY", "ZIIXYYYIII", "YZYXZZYYYI", "XIXXXXZIYI", "IIIYXIZYYI", 
        "YZIZIYXYII", "XXYXYXYXII", "YIYZYIZYXY", "XZIZIIIIYI", "XXIXIZZXXY", "IZIIYYZIXY", 
        "ZZZZYIXYZY", "IXXXZXZZXY", "ZZYIXXIXII", "YYXYIIIZII", "YIIIXZIZII", "YIZXYZXIXY", 
        "ZIIXIXIXZY", "IXZYYXZXYI", "ZZYXYIYYYI", "YYXZXYZYXY", "IIIIZXYYYI", "YYIYIXYIZY", 
        "YIXIZIZYZY", "YYIIXYXYXY", "YZIZXIIIII", "XYXIXYYZII", "YIXZXXIYYI", "XZZZYIXYXY", 
        "IXYYIIIIYI", "YXYIYIYIZY", "IXYIYXXIYI", "XXYIZXYIYI", "IXZIYXXYZY", "ZIYIZIIYXY", 
        "XZZIXIXXXY"
    ]
    starterList=["YZZYZXYZII"]
    
    extra=[]
    ~,L= GenerateGroup!(starterList, QubitNumber)
    Len=[L]

    for i=2:length(starters)
        push!(starterList,starters[i])
        ~,L= GenerateGroup!(starterList, QubitNumber)
        @show i
        if Len[end] == L
            push!(extra,starters[i])
            deleteat!(starterList, findall(x->x==starters[i],starterList))
        else
            push!(Len,L)
        end
    end

    Group,L= GenerateGroup!(starterList, QubitNumber)
    for i=1:length(Group)
        @show CheckInseprability([Group[i];extra[1]])
    end

    
    QubitNumber=2
    Starters=["I","X","Y","Z"]
    Starters=["IX","IY","IZ","II","XI","XY","XZ","XX","YI","YX","YY","YZ","ZI","ZX","ZY","ZZ"]
    GeneratePauliStrings(Starters)

end
