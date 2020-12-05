import os

def get_mess(file_name, begin_time, end_time):
    file = open(file_name, "r+")
    all_line = file.readlines()
    message = ""
    now = begin_time
    has_head = True if begin_time==1 else False

    for i in range(len(all_line)):
        words = all_line[i].split(" ")
        #print(words[0])
        if words[0] == "MODEL":
            message += "MODEL       " + str(now) + "\n"
            now += 1
            has_head = True
            #print("nnn")
        elif words[0] == "ENDMDL" or words[0] == "ENDMDL\n":
            message += "ENDMDL\n"
            print(now)
            if now >= end_time:
                break
        else:
            if not has_head:
                continue
            message += all_line[i]
    #print(message)
    return message
def get_pdb(pdb_numb, rank, pdb_iteration, n_iteration):
    #pdb_iteration 存pdb步长
    #n_iteration 总步长
    pdb_file_name = ["out" + str(rank) + "-" + str(id) + ".pdb" for id in range(pdb_numb)]
    pdb_file_name[0] = "out" + str(rank) + ".pdb"
    #print(pdb_file_name)
    restart_file = open("restart.log",'r+')
    merge_pdb = open("all%i.pdb"%(rank), 'w')
    #get every pdb begin time
    pdb_begin_time = [id for id in range(pdb_numb)]
    all_line = restart_file.readlines()
    for id in range(len(all_line)):
        pdb_begin_time[id] = int(all_line[id])
    pdb_begin_time[pdb_numb - 1] = n_iteration
    print(pdb_begin_time)
    begin_time = 0
    for id in range(pdb_numb):
        message = get_mess(pdb_file_name[id], begin_time * pdb_iteration + 1 , (pdb_begin_time[id] * pdb_iteration + 1))
        begin_time = pdb_begin_time[id]
        merge_pdb.write(message)

n_replicas = 4
n_iterations = 100
savelogstep = 500
savepdbstep = 500
ex_interval = 1000
pdb_numb = 2;
pdb_iteration = ex_interval/savepdbstep
n_iteration = n_iterations
for i in range(n_replicas):
    get_pdb(pdb_numb,i,pdb_iteration, n_iteration)



